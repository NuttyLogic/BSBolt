from BSB.BSB_Align.AlignmentHelpers import convert_alpha_numeric_cigar, get_mapping_length, convert_cigar_tuple
from BSB.BSB_Utils.UtilityFunctions import reverse_complement


class ProcessSamAlignment:
    """

    """

    def __init__(self, contig_lens=None):
        self.contig_lens = contig_lens
        # when crick strands are reversed, mapping strand flips so set appropriate sam flag
        self.sense_flags = {'0', '99', '163', '256', '323', '419', '65', '67', '73', '97', '99', '129', '137', '161',
                            '163', '321', '323', '329', '353', '355', '385', '393', '417', '419'}
        self.mixed_unmapped = {'69', '101', '133', '165'}
        self.watson_mapping_conversion = {'0': '0', '16': '0', '256': '256', '272': '256',
                                          '99': '67', '147': '131', '83': '67', '163': '131',
                                          '355': '323', '403': '387', '339': '323', '419': '387',
                                          '97': '65', '145': '129', '81': '65', '161': '129',
                                          '65': '65', '129': '129', '113': '65', '177': '129',
                                          '353': '321', '401': '385', '337': '321', '417': '385',
                                          '321': '321', '385': '385', '369': '321', '433': '385',
                                          '73': '73', '329': '329', '133': '133',
                                          '137': '137', '393': '393', '69': '69',
                                          '89': '73', '345': '329', '165': '133',
                                          '153': '137', '409': '393', '101': '69'}
        self.crick_mapping_conversion = {'0': '16', '16': '16', '256': '272', '272': '272',
                                         '99': '115', '147': '179', '83': '115', '163': '179',
                                         '355': '371', '403': '435', '339': '371', '419': '435',
                                         '97': '113', '145': '177', '81': '113', '161': '177',
                                         '65': '113', '129': '177', '113': '113', '177': '177',
                                         '353': '369', '401': '433', '337': '369', '417': '433',
                                         '321': '369', '385': '433', '369': '369', '433': '433',
                                         '73': '89', '329': '345', '133': '165',
                                         '137': '153', '393': '409', '69': '101',
                                         '89': '89', '345': '345', '165': '165',
                                         '153': '153', '409': '409', '101': '101'}

    def process_read(self, sam_read):
        sense = sam_read['FLAG'] in self.sense_flags
        mixed_unmapped = sam_read['FLAG'] in self.mixed_unmapped
        if 'C_' in sam_read['mapping_reference']:
            self.format_crick_reads(sam_read, sense, mixed_unmapped)
        else:
            self.format_watson_reads(sam_read, sense, mixed_unmapped)
        assert len(sam_read['SEQ']) == len(sam_read['QUAL']), "Length of read sequence != Length of read " \
                                                              "quality, check file formatting "
        sam_read['SAM_TAGS'].append(f'XO:Z:{sam_read["mapping_reference"]}')

    def format_watson_reads(self, sam_read, sense, mixed_unmapped):
        """Format Watson reads so all reads are on positive strand, change flags to represent this"""
        if not sense and not mixed_unmapped:
            sam_read['SEQ'] = reverse_complement(sam_read['SEQ'])
            sam_read['QUAL'] = sam_read['QUAL'][::-1]
        sam_read['FLAG'] = self.watson_mapping_conversion[sam_read['FLAG']]

    def format_crick_reads(self, sam_read, sense, mixed_unmapped):
        """Format crick reads, some reads reverse complemented others only reversed"""
        # if sequence is mapped to -crick, reverse complement sequence
        if sense and not mixed_unmapped:
            sam_read['SEQ'] = reverse_complement(sam_read['SEQ'])
            sam_read['QUAL'] = sam_read['QUAL'][::-1]
        # reverse cigar string
        cigar_tuple = convert_alpha_numeric_cigar(sam_read['CIGAR'])[::-1]
        # convert flag to be relative to negative strand
        sam_read['FLAG'] = self.crick_mapping_conversion[sam_read['FLAG']]
        # get mapping length of the read to adjust coordinates relative to watson strand
        mapping_length = get_mapping_length(cigar_tuple)
        mapping_loc = self.contig_lens[sam_read['RNAME']] - int(sam_read['POS']) - mapping_length + 2
        sam_read['POS'] = str(mapping_loc)
        if sam_read['PNEXT'] != '*':
            mate_chrom = sam_read['RNAME'] if sam_read['RNEXT'] == '=' else sam_read['RNEXT']
            mate_pos = self.contig_lens[mate_chrom] - int(sam_read['PNEXT']) - mapping_length + 2
            sam_read['PNEXT'] = str(mate_pos)
        # convert cigar back to alpha numeric representation for mapped reads
        if not mixed_unmapped:
            sam_read['CIGAR'] = convert_cigar_tuple(cigar_tuple)
        # change template length for paired end reads
        if sam_read['TLEN'] != '0':
            tlen = sam_read['TLEN']
            sam_read['TLEN'] = f'-{tlen}' if '-' not in tlen else tlen.replace('-', '')
