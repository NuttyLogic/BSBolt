import pickle
from typing import Any, Dict, List, Set
import pysam
from BSBolt.Align.LaunchBowtie2Alignment import Bowtie2Align
from BSBolt.Align.ProcessSamReads import ProcessSamAlignment
from BSBolt.Align.AlignmentHelpers import write_bam_line
from BSBolt.Align.SortReads import SortReads


class BisulfiteAlignmentAndProcessing:
    """ Unifying alignment and processing class. Class launches bowtie2 alignment. Alignment call is blocking, so class
    pauses until external alignment until complete. If there is an alignment error, error is propagated through the
    class. Following alignment, fastq reads and temporary sam lines, converted to tab format, are iterated through
    together, processed, and written to final sam file. Paired end and single end specific processing happens within
    class, individual read processing is handled by ProcessSamAlignment class.
    Keyword Arguments:
        fastq1 (str): path to fastq file
        fastq2 (str): path to fastq1 mate pair
        undirectional_library (bool): perform undirectional alignment
        bowtie2_commands (list of str): commands to use for bowtie2 alignment
        bsb_database(str): path to bsb 2 directory
        bowtie2_path(str): path to bowtie2 executable
        output_path(str): output prefix of output file
        mismatch_threshold (int): mismatch threshold to consider a valid read
        command_line_arg (str): list of commandline args to output with sam file
        non_converted_output (bool): map non in silico converted reads to bisulfite reference
        allow_discordant (bool): all discordant and mixed read pairs
    Attributes:
        self.bowtie2_mapping (dict): argument for launching bowtie2 mapping
        self.undirectional_library (bool): perform undirectional alignment
        self.bsb_database(str): path to bsb database directory
        self.output_path (str): output prefix of output file
        self.mismatch_threshold (int): mismatch threshold to consider a valid read
        self.command_line_arg (str): list of commandline args to output with sam file
        self.sam_output (pysam.AlignmentFile): TextIO instance to output processed sam reads.
        self.mapping_statistics (dict): dict of mapping statistics
        self.read_sorter (SortReads):  sort and return mapped / unmapped reads
        self.read_processor (ProcessSamAlignment): process mapped reads and orient them to reference strand
    """

    def __init__(self, fastq1: str = None, fastq2: str = None, undirectional_library=False,
                 bowtie2_commands: List[str] = None, bsb_database: str = None, bowtie2_path: str = None,
                 output_path: str = None, mismatch_threshold=4, command_line_arg: str = None,
                 non_converted_output=False, allow_discordant=False):
        self.bowtie2_mapping = dict(fastq1=fastq1, fastq2=fastq2, undirectional_library=undirectional_library,
                                    bowtie2_commands=bowtie2_commands, bowtie2_path=bowtie2_path,
                                    bsb_database=bsb_database, no_conversion=non_converted_output)
        self.bsb_database = bsb_database
        self.command_line_arg = command_line_arg
        self.output_path = output_path
        self.contig_lens = self.get_contig_lens
        self.sam_output = self.get_output_object
        self.mapping_statistics = dict(total_reads=0,
                                       unmapped_reads=0, multireference_reads=0,
                                       reads_mapped_1=0, reads_mapped_more_than_1=0,
                                       discordant_reads_1=0, discordant_reads_more_than_1=0,
                                       mixed_reads_1=0, mixed_reads_more_than_1=0,
                                       W_C2T=0, W_G2A=0,
                                       C_C2T=0, C_G2A=0)
        self.read_sorter = SortReads(mismatch_threshold=mismatch_threshold,
                                     allow_discordant=allow_discordant,
                                     contig_lens=self.contig_lens)
        self.read_processor = ProcessSamAlignment()

    @property
    def get_output_object(self) -> pysam.AlignmentFile:
        """
        Returns: TextIO instance to output processed sam reads
        """
        # load contig len dictionary
        with open(f'{self.bsb_database}genome_index.pkl', 'rb') as index:
            genome_index: dict = pickle.load(index)
        sam_header = {'HD': {'VN': '1.0'}, 'SQ': []}
        # add chromosome information to bam header
        for chromosome, chromosome_length in genome_index.items():
            sam_header['SQ'].append({'LN': chromosome_length, 'SN': chromosome})
        sam_header['PG'] = [{'ID': '01', 'PN': 'BSBolt', 'CL': self.command_line_arg, 'VN': '0.0.2'}]
        sam_out = pysam.AlignmentFile(f'{self.output_path}.bam', 'wb', header=sam_header)
        return sam_out

    def align_reads(self):
        """ Instance Bowtie2Alingment class and map reads"""
        alignment = iter(Bowtie2Align(**self.bowtie2_mapping))
        sam_read = next(alignment)
        read_group, read_name = [sam_read], sam_read['QNAME']
        # iterate and process reads
        while True:
            try:
                sam_read = next(alignment)
            except StopIteration:
                self.process_reads(read_group)
                break
            else:
                if sam_read['QNAME'] != read_name:
                    self.process_reads(read_group)
                    read_group, read_name = [sam_read], sam_read['QNAME']
                else:
                    read_group.append(sam_read)

    def process_reads(self, read_group: List[Dict[str, Any]]):
        # retrieve sorted reads
        read_status, mapping_references, sorted_reads = self.read_sorter.get_reads(read_group=read_group)
        # processes sorted reads and update number of observed reads by read name
        self.process_sam_reads(read_status, mapping_references, sorted_reads)
        self.mapping_statistics['total_reads'] += 1

    @property
    def get_contig_lens(self) -> Dict[str, int]:
        with open(f'{self.bsb_database}genome_index.pkl', 'rb') as contig_lens:
            return pickle.load(contig_lens)

    def update_mapping_statistics(self, read_status, mapping_references, sorted_reads):
        """Mapping statistics are updated per read group based on the mapping status"""
        if not read_status:
            if mapping_references:
                self.mapping_statistics['multireference_reads'] += 1
            else:
                self.mapping_statistics['unmapped_reads'] += 1
        elif read_status in {'Single', 'Paired'}:
            self.mapping_statistics[sorted_reads[0]['mapping_reference']] += 1
            read_len = len(sorted_reads)
            if read_status == 'Paired':
                read_len = read_len / 2
            if read_len > 1:
                self.mapping_statistics['reads_mapped_more_than_1'] += 1
            else:
                self.mapping_statistics['reads_mapped_1'] += 1
        elif read_status == 'Discordant':
            self.mapping_statistics[sorted_reads[0]['mapping_reference']] += 1
            if len(sorted_reads) > 2:
                self.mapping_statistics['discordant_reads_more_than_1'] += 1
            else:
                self.mapping_statistics['discordant_reads_1'] += 1
        else:
            self.mapping_statistics[sorted_reads[0]['mapping_reference']] += 1
            if len(sorted_reads) > 2:
                self.mapping_statistics['mixed_reads_more_than_1'] += 1
            else:
                self.mapping_statistics['mixed_reads_1'] += 1

    def process_sam_reads(self, read_status: bool, mapping_references: Set[str],
                          sorted_reads: List[Dict[str, Any]]):
        """Launch sam read processing
        Arguments:
            read_status (str): type of reads being passed, (single, paired, mapped, discordant, mixed)
            mapping_references (list or None): list of mapping strands
            sorted_reads (List[Dict[str, Union[str, List, int]]]): list of processed sam reads"""
        self.update_mapping_statistics(read_status, mapping_references, sorted_reads)
        for read in sorted_reads:
            if read_status:
                self.read_processor.process_read(sam_read=read)
            write_bam_line(read, self.sam_output)
