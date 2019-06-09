import pickle
import pysam
from BSB.BSB_Align.LaunchBowtie2Alignment import Bowtie2Align
from BSB.BSB_Align.ProcessSamReads import ProcessSamAlignment
from BSB.BSB_Align.AlignmentHelpers import write_bam_line


class BisulfiteAlignmentAndProcessing:
    """ Unifying alignment and processing class. Class launches bowtie2 alignment. Alignment call is blocking, so class
    pauses until external alignment until complete. If there is an alignment error, error is propagated through the
    class. Following alignment, fastq reads and temporary sam lines, converted to tab format, are iterated through
    together, processed, and written to final sam file. Paired end and single end specific processing happens within
    class, individual read processing is handled by ProccessSamAlingment class.
    Keyword Arguments:
        fastq1 (str): path to fastq file
        fastq2 (str): path to fastq1 mate pair
        undirectional_library (bool): perform undirectional alignment
        bowtie2_commands (list of str): commands to use for bowtie2 alignment
        bsb_database(str): path to bsb 2 directory
        bowtie2_path(str): path to bowtie2 executable
        output_path(str): output prefix of output file
        conversion_threshold (tuple (float, int)): proportion of unconverted CHH sites and minimum number of CHH sites
                                                    to label a read incompletely converted
        mismatch_threshold (int): mismatch threshold to consider a valid read
        command_line_arg (str): list of commandline args to output with sam file
    Attributes:
        self.bowtie2_mapping (dict): argument for launching bowtie2 mapping
        self.fastq1 (str): path to fastq file
        self.fastq2 (str): path to fastq1 mate pair
        self.undirectional_library (bool): perform undirectional alignment
        self.bsb_database(str): path to bsb database directory
        self.paired_end (bool): paired end reads
        self.output_path (str): output prefix of output file
        self.self.sam_tuple (tuple): tuple of temporary output files
        self.conversion_threshold (tuple (float, int)): proportion of unconverted CHH sites and minimum number of
                                                        CHH sites to label a read incompletely converted
        self.mismatch_threshold (int): mismatch threshold to consider a valid read
        self.command_line_arg (str): list of commandline args to output with sam file
        self.contig_sequence_dict (dict): dict of lengths, bp, of contig sequences, populated in ProccessSamAlingment
                                          class instances
        self.sam_output (TextIO): TextIO instance to output processed sam reads.
        self.mapping_statistics (dict): dict of mapping statistics
        self.flag_correction (dict): dict containing specific temp sam file processing instructions
    """

    def __init__(self, fastq1=None, fastq2=None, undirectional_library=False, bowtie2_commands=None,
                 bsb_database=None, bowtie2_path=None, output_path=None,
                 mismatch_threshold=4, command_line_arg=None, unmapped_output=False):
        self.bowtie2_mapping = dict(fastq1=fastq1, fastq2=fastq2, undirectional_library=undirectional_library,
                                    bowtie2_commands=bowtie2_commands, bowtie2_path=bowtie2_path,
                                    bsb_database=bsb_database)
        self.bsb_database = bsb_database
        self.undirectional_library = undirectional_library
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.command_line_arg = command_line_arg
        self.unmapped_output = unmapped_output
        self.paired_end = False
        # if self.fastq2 provided assume library is paired end
        if self.fastq2:
            self.paired_end = True
        self.output_path = output_path
        self.mismatch_threshold = mismatch_threshold
        self.contig_lens = self.get_contig_lens
        self.sam_output = self.get_output_object
        self.mapping_statistics = dict(total_reads=0, multimapped_reads=0,
                                       unmapped_reads=0, multireference_reads=0,
                                       unique_reads=0,
                                       W_C2T=0, W_G2A=0,
                                       C_C2T=0, C_G2A=0)

    @property
    def get_output_object(self):
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
        sam_reads, read_name = [], None
        while True:
            self.mapping_statistics['total_reads'] += 1
            try:
                sam_read = next(alignment)
            except StopIteration:
                break
            if not read_name:
                read_name = sam_read['QNAME']
            elif sam_read['QNAME'] != read_name:
                self.process_sam_reads(sam_reads)
                sam_reads, read_name = [], sam_read['QNAME']
            if self.paired_end:
                sam_read_paired = next(alignment)
                assert sam_read['QNAME'] == sam_read_paired['QNAME'], '.fastq files are not paired, sort .fastq files'
                sam_reads.append((sam_read, sam_read_paired))
            else:
                sam_reads.append(tuple(sam_read))
        self.process_sam_reads(sam_reads)

    @property
    def get_contig_lens(self):
        with open(f'{self.bsb_database}genome_index.pkl', 'rb') as contig_lens:
            return pickle.load(contig_lens)

    def write_alignment_reads(self, read_grouping, output_object):
        for read in read_grouping:
            write_bam_line(read, self.sam_output)

    def process_sam_reads(self, sam_reads):
        """Launch sam read processing
        Arguments:
            sam_reads (list): list of sam_read dictionaries and fastq lists"""
        process_sam = ProcessSamAlignment(sam_reads=sam_reads,
                                          contig_lens=self.contig_lens,
                                          mismatch_threshold=self.mismatch_threshold)
        mapping_number, processed_reads = process_sam.output_reads
        if mapping_number == 1:
            for read_grouping in processed_reads:
                self.write_alignment_reads(read_grouping, self.sam_output)
            self.mapping_statistics[processed_reads[0][1]['mapping_reference']] += 1
            if len(processed_reads) > 1:
                self.mapping_statistics['multimapped_reads'] += 1
            else:
                self.mapping_statistics['unique_reads'] += 1
        else:
            if self.unmapped_output:
                for read_grouping in processed_reads:
                    self.write_alignment_reads(read_grouping, self.sam_output)
            if mapping_number < 1:
                self.mapping_statistics['unmapped_reads'] += 1
            else:
                self.mapping_statistics['multireference_reads'] += 1

