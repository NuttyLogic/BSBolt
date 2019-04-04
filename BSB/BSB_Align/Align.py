import pickle
import subprocess
import pysam
from BSB.BSB_Align.BT2_Mapping import Bowtie2Alignment
from BSB.BSB_Align.ProcessSamReads import ProcessSamAlignment
from BSB.BSB_Align.TabSamIterator import TabSamIterator
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
                 bsb_database=None, bowtie2_path=None, output_path=None, conversion_threshold=(0.5, 5),
                 mismatch_threshold=None, command_line_arg=None, non_converted_output=False):
        self.bowtie2_mapping = dict(fastq1=fastq1, fastq2=fastq2, undirectional_library=undirectional_library,
                                    bowtie2_commands=bowtie2_commands, bowtie2_path=bowtie2_path,
                                    output_path=output_path, bsb_database=bsb_database,
                                    non_converted_output=non_converted_output)
        assert isinstance(undirectional_library, bool)
        self.bsb_database = bsb_database
        self.undirectional_library = undirectional_library
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.command_line_arg = command_line_arg
        self.paired_end = False
        # if self.fastq2 provided assume library is paired end
        if self.fastq2:
            self.paired_end = True
        self.output_path = output_path
        self.sam_tuple: tuple = self.get_sam_tuple
        self.conversion_threshold = conversion_threshold
        self.mismatch_threshold = mismatch_threshold
        self.contig_sequence_dict = {}
        self.sam_output = self.get_output_object
        self.mapping_statistics = dict(total_reads=0, multimapped_reads=0, unmapped_reads=0)
        self.flag_correction = {'W_C2T': {'99': '67', '147': '131'},
                                'W_G2A': {'99': '115', '147': '179'},
                                'C_C2T': {'99': '179', '147': '115'},
                                'C_G2A': {'99': '131', '147': '67'}}

    @property
    def get_sam_tuple(self):
        # alignment for one strand pair if library is undirectional, all if directional
        sam_list = [f'{self.output_path}.W_C2T.sam.temp', f'{self.output_path}.C_C2T.sam.temp']
        if self.undirectional_library:
            sam_list.extend([f'{self.output_path}.W_G2A.sam.temp', f'{self.output_path}.C_G2A.sam.temp'])
        return tuple(sam_list)

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

    def launch_bisulfite_aligment(self):
        """ Instance Bowtie2Alingment class and map reads"""
        bisulfite_alignment = Bowtie2Alignment(**self.bowtie2_mapping)
        bisulfite_alignment.launch_bt2_mapping()

    def process_reads(self):
        """ Iterate through fastq files and temp sam files generated with bowtie2 mapping and combine to one
        sam file"""
        # instance iterator to go through fastq lines and sam lines together, this assumes the order is identical for
        # all alignment files (if this isn't true alignment processing will quickly break)
        tab_sam_iterator = iter(TabSamIterator(fastq1=self.fastq1, fastq2=self.fastq2, sam_tuple=self.sam_tuple))
        while True:
            # get next line or break
            try:
                sam_reads = next(tab_sam_iterator)
            except StopIteration:
                break
            else:
                # get sam read information for first reads
                mapping_number, processed_read, bisulfite_strand = self.process_sam_reads(sam_reads)
                # if paired end get second read information and perform PE specific processing
                if self.paired_end:
                    # if pe sam_read_2 should always be present
                    sam_reads_2 = next(tab_sam_iterator)
                    mapping_number_2, processed_read_2, bisulfite_strand_2 = self.process_sam_reads(sam_reads_2)
                    # only perform additional processing if both reads uniquely aligned and passed filters
                    if bisulfite_strand_2 and bisulfite_strand:
                        # sanity check to make sure reads are paired end
                        assert processed_read['QNAME'].split('/')[0] == processed_read_2['QNAME'].split('/')[0]
                        # depending on conversion and processing correction make location may be incorrect so update
                        processed_read['PNEXT'] = str(processed_read_2['POS'])
                        processed_read_2['PNEXT'] = str(processed_read['POS'])
                        # set proper template sign if position updated during read processing
                        self.check_template(processed_read, processed_read_2)
                        # set proper sam flags
                        self.correct_paired_end_flag(bisulfite_strand, processed_read, processed_read_2)
                        # output lines
                        self.output_sam_lines(bisulfite_strand, processed_read, bisulfite_strand_2, processed_read_2)
                    self.update_mapping_statistics(mapping_number_2, bisulfite_strand_2)
                else:
                    # output se line
                    self.output_sam_lines(bisulfite_strand, processed_read)
            self.update_mapping_statistics(mapping_number, bisulfite_strand)
        # close output file
        self.sam_output.close()

    def output_sam_lines(self, bisulfite_strand, processed_read, bisulfite_strand_2=None, processed_read_2=None):
        """ Write out sam files if read mapped uniquely and pass quality filters. Handles PE and SE reads.
        Arguments
            bisulfite_strand (str or None): mapping strand
            processed_read (dict): read dictionary
        Keyword Arguments:
            bisulfite_strand_2 (str or None): mapping strand
            processed_read_2 (dict): read dictionary
        """
        # if process_read_2 given output is PE
        if processed_read_2:
            # both reads have to map uniquely to write out to sam file
            if bisulfite_strand and bisulfite_strand_2:
                write_bam_line(processed_read, self.sam_output)
                write_bam_line(processed_read_2, self.sam_output)
        else:
            # single end output
            if bisulfite_strand:
                # correct single end flag
                self.correct_single_end_flag(bisulfite_strand, processed_read)
                write_bam_line(processed_read, self.sam_output)

    def clean_temp_files(self):
        for sam_file in self.sam_tuple:
            subprocess.run(['rm', sam_file])

    @staticmethod
    def correct_single_end_flag(bisulfite_strand, processed_read):
        """Change flag for reverse complement strands"""
        if bisulfite_strand == 'W_G2A' or bisulfite_strand == 'C_C2T':
            processed_read['FLAG'] = '16'

    def process_sam_reads(self, sam_reads):
        """Launch sam read processing
        Arguments:
            sam_reads (list): list of sam_read dictionaries and fastq lists"""
        process_sam = ProcessSamAlignment(sam_line_dict=sam_reads,
                                          contig_dict=self.contig_sequence_dict,
                                          bsb_database=self.bsb_database,
                                          conversion_threshold=self.conversion_threshold,
                                          mismatch_threshold=self.mismatch_threshold)
        mapping_number, processed_read, bisulfite_strand = process_sam.output_read
        # return results
        return mapping_number, processed_read, bisulfite_strand

    def correct_paired_end_flag(self, bisulfite_strand, processed_read, processed_read_2):
        """Properly orient flags based on mapping strand
        Arguments:
            bisulfite_strand (str or None): mapping strand
            processed_read (dict): read dictionary
            processed_read_2 (dict): read dictionary """
        processed_read['FLAG'] = self.flag_correction[bisulfite_strand][processed_read['FLAG']]
        processed_read_2['FLAG'] = self.flag_correction[bisulfite_strand][processed_read_2['FLAG']]

    @staticmethod
    def check_template(read_1, read_2):
        """Set proper template lengths, if upstream and downstream reads are flipped during processing template length
        signs all need to be flipped
        Arguments:
            read_1 (dict): read dictionary
            read_2 (dict): read dictionary """
        # check read position and if read position flipped flip template length
        if int(read_1['POS']) > int(read_2['POS']):
            if int(read_1['TLEN']) > 0:
                tlen_1 = str(read_1['TLEN'])
                read_1['TLEN'] = str(read_2['TLEN'])
                read_2['TLEN'] = tlen_1
        else:
            # if template length sign incorrect flip sign
            if int(read_1['TLEN']) < 0:
                tlen_1 = str(read_1['TLEN'])
                read_1['TLEN'] = str(read_2['TLEN'])
                read_2['TLEN'] = tlen_1

    def update_mapping_statistics(self, mapping_number, bisulfite_strand):
        """ Update mapping statistics based on the number of time read mapped
        Arguments:
            mapping_number (int): number of times read mapped
            bisulfite_strand (str or none): if one mapping bisulfite strand else None"""
        self.mapping_statistics['total_reads'] += 1
        if bisulfite_strand:
            try:
                self.mapping_statistics[bisulfite_strand] += 1
            except KeyError:
                self.mapping_statistics[bisulfite_strand] = 1
        else:
            if mapping_number < 1:
                self.mapping_statistics['unmapped_reads'] += 1
            else:
                self.mapping_statistics['multimapped_reads'] += 1
