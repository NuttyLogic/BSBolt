from BSB.BSB_Align.AlignmentHelpers import convert_alpha_numeric_cigar, get_mapping_length, convert_cigar_tuple
from BSB.BSB_Utils.UtilityFunctions import reverse_complement


class ProcessSamAlignment:
    """ SAM alignment processing class. Takes combined fastq and sam reads. Outputs process SAM reads, if read is
    uniquely aligned the correctly mapped read will be sorted. If read is multimapped or unmapped the original fastq
    sequence will be output.
    Keyword Arguments:
        sam_reads (list): list of paired/unpaired reads from bowtie2 alignment
        mismatch_threshold (int): number of mismatches allowed for a read to be considered
        contig_lens (dict): lens of mapping contig to reverse crick mapping reads
    Attributes:
        self.sam_reads (list): list of reads paired/unpaired from bowtie2
        self.contig_lens (list): lens of mapping contigs to reverse crick mapping
        self.mismatch_threshold (int): number of mismatches allowed for a read to be considered
        self.mapping_flags (set): sam flags that indicate a properly aligned read or read pair
        self.output_read (object): processed sam read handled by BSB_Align
        self.watson_mapping_conversion (dict_: appropriate flags after changing watson strand
        self.crick_mapping_conversion (dict): appropriate flags after changing crick strand  
    """

    def __init__(self, sam_reads=None, contig_lens=None, mismatch_threshold=4):
        self.sam_reads = sam_reads
        self.mismatch_threshold = mismatch_threshold
        self.contig_lens = contig_lens
        self.mapping_flags = {'0', '16', '83', '99', '147', '163', '256', '272',
                              '339', '355', '403', '419'}
        # when crick strands are reversed, mapping strand flips so set appropriate sam flag
        self.watson_mapping_conversion = {'0': '0', '16': '0', '256': '256', '272': '256',
                                          '99': '67', '147': '131', '83': '67', '163': '131',
                                          '355': '323', '403': '387', '339': '323', '419': '387'}
        self.crick_mapping_conversion = {'0': '16', '16': '16', '256': '272', '272': '272',
                                         '99': '115', '147': '179', '83': '115', '163': '179',
                                         '355': '371', '403': '435', '339': '371', '419': '435'}
        self.output_reads = self.get_combined_bisulfite_read

    @property
    def get_combined_bisulfite_read(self):
        # get number of correctly maped reads
        mapped_read_number, mapped_reads = self.check_read_mapping
        if mapped_read_number == 1:
            processed_reads = []
            for read_grouping in mapped_reads:
                reads = [self.process_read(read, len(read_grouping)) for read in read_grouping]
                if len(reads) > 1:
                    read_1_pos, read_2_pos = str(reads[0]['POS']), str(reads[1]['POS'])
                    reads[0]['PNEXT'], reads[1]['PNEXT'] = read_2_pos, read_1_pos
                processed_reads.append(tuple(reads))
            return mapped_read_number, processed_reads
        return mapped_read_number, mapped_reads

    @property
    def check_read_mapping(self):
        # processing assumes standard input order
        mapped_strands = set()
        mapping_reads = []
        for read_group in self.sam_reads:
            if read_group[0]['FLAG'] in self.mapping_flags:
                mismatch_check = False
                for read in read_group:
                    if self.check_mismatch(read) or read['FLAG'] not in self.mapping_flags:
                        mismatch_check = True
                if not mismatch_check:
                    mapped_strands.add(read_group[0]['mapping_reference'])
                    mapping_reads.append(read_group)
        if len(mapped_strands) == 1:
            return len(mapped_strands), mapping_reads
        return len(mapped_strands), [self.set_unmapped_read(self.sam_reads[0], mapped_strands)]

    @staticmethod
    def set_unmapped_read(read_group, mapping_strands):
        unmapped_flags = tuple('4') if len(read_group) == 1 else ('77', '141')
        mapping_strand_id = ','.join(mapping_strands)
        unmapped_reads = []
        for read, flag in zip(read_group, unmapped_flags):
            unmapped_read = dict(read)
            unmapped_read['FLAG'] = flag
            unmapped_read['RNAME'] = '*'
            unmapped_read['POS'] = '0'
            unmapped_read['CIGAR'] = '*'
            unmapped_read['PNEXT'] = '0'
            unmapped_read['TLEN'] = '0'
            if mapping_strand_id:
                unmapped_read['SAM_TAGS'] = ['YT:Z:UP', f'XO:Z:{mapping_strand_id}']
            else:
                unmapped_read['SAM_TAGS'] = ['YT:Z:UP']
            unmapped_reads.append(unmapped_read)
        return tuple(unmapped_reads)

    def check_mismatch(self, read):
        """ Arguments:
                read (dict): individual sam read
        """
        # retrieve sam tags
        SAM_TAGS: list = read['SAM_TAGS']
        for tag in SAM_TAGS:
            # if sam tag XM, mismatch sam tag
            if 'XM:i' in tag:
                # parse sam tag
                mismatch_number = int(tag.split(':')[-1])
                if mismatch_number <= self.mismatch_threshold:
                    return False
                return True
        return True

    def process_read(self, read, paired_end):
        mapped_read = dict(read)
        contig_len = self.contig_lens[mapped_read['RNAME']]
        # rule to reverse flag based on mapping strand
        if 'C_' in mapped_read['mapping_reference']:
            self.format_crick_reads(mapped_read, contig_len, paired_end)
        else:
            self.format_watson_reads(mapped_read)
        assert len(mapped_read['SEQ']) == len(mapped_read['QUAL']), "Length of read sequence != Length of read " \
                                                                    "quality, check file formatting "
        mapped_read['SAM_TAGS'].append(f'XO:Z:{mapped_read["mapping_reference"]}')
        return mapped_read

    def format_watson_reads(self, mapped_read):
        """Format Watson reads so all reads are on positive strand, change flags to represent this"""
        complement_flags = {'16', '83', '147', '272', '339', '403'}
        if mapped_read['FLAG'] in complement_flags:
            mapped_read['SEQ'] = reverse_complement(mapped_read['SEQ'])
            cigar_tuple = convert_alpha_numeric_cigar(mapped_read['CIGAR'])
            mapped_read['CIGAR'] = convert_cigar_tuple(cigar_tuple[::-1])
        mapped_read['FLAG'] = self.watson_mapping_conversion[mapped_read['FLAG']]

    def format_crick_reads(self, mapped_read, contig_len, paired_end):
        """Format crick reads, some reads reverse complemented others only reversed"""
        reverse_complement_flags = {'0', '99', '163', '256', '323', '419'}
        # cigar tuple used for both strand
        cigar_tuple = convert_alpha_numeric_cigar(mapped_read['CIGAR'])
        # if sequence is mapped to -crick, reverse complement sequence otherwise just reverse the sequence
        if mapped_read['FLAG'] in reverse_complement_flags:
            mapped_read['SEQ'] = reverse_complement(mapped_read['SEQ'])
            cigar_tuple = convert_alpha_numeric_cigar(mapped_read['CIGAR'])[::-1]
            mapped_read['QUAL'] = mapped_read['QUAL'][::-1]
        # convert flag to be relative to negative strand
        mapped_read['FLAG'] = self.crick_mapping_conversion[mapped_read['FLAG']]
        # get mapping length of the read to adjust coordinates relative to watson strand
        mapping_length = get_mapping_length(cigar_tuple)
        mapping_loc = contig_len - int(mapped_read['POS']) - mapping_length + 2
        mapped_read['POS'] = str(mapping_loc)
        # convert cigar back to alpha numeric representation
        mapped_read['CIGAR'] = convert_cigar_tuple(cigar_tuple)
        # change template length for paired end reads
        if paired_end > 1:
            tlen = mapped_read['TLEN']
            mapped_read['TLEN'] = f'-{tlen}' if '-' not in tlen else tlen.replace('-', '')
