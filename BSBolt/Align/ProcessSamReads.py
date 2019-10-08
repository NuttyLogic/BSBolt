from BSBolt.Utils.UtilityFunctions import reverse_complement


class ProcessSamAlignment:
    """
    Sets mapped reads relative to the regular genomic space by placing all watson reads on sense strand and placing
    all crick reads on anti-sense strand. As crick reads map to reverse complement strand coordinate space has to
    adjusted. Note, unmapped reads with a mapped mate aren't modified except for the report flag to reflect the
    mapping strand of the mate.
    Keyword Arguments:
        self.contig_lens (dict[str: int]): chromosome mapping to the length of the mapping contig used to
                                           put crick reads into normal mapping space
    Attributes:
        self.sense_flags (set): flags associated with reads mapped to the sense strand
        self.mixed_unmapped (set): flags for unmapped reads that don't have to be modified unless the
                                   mate read switches strands
        self.watson_mapping_conversion (dict): flag changes to move reads relative to sense strand
        self.crick_mapping_conversion (dict): flag changes to move reads relative to the anti-sense strand
    """

    def __init__(self):
        # when crick strands are reversed, mapping strand flips so set appropriate sam flag
        self.sense_flags = {'0', '99', '163', '256', '323', '419', '65', '67', '73', '97', '99', '129', '131',
                            '137', '161', '163', '321', '323', '329', '353', '355', '385', '393', '417', '419'}
        self.mixed_unmapped = {'69', '101', '133', '165'}
        self.watson_mapping_conversion = {'0': '0', '16': '0', '256': '256', '272': '256',
                                          '99': '67', '147': '131', '83': '67', '163': '131',
                                          '67': '67', '131': '131', '115': '67', '179': '131',
                                          '355': '323', '403': '387', '339': '323', '419': '387',
                                          '323': '323', '387': '387', '371': '323', '435': '387',
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
                                         '67': '115', '131': '179', '115': '115', '179': '179',
                                         '355': '371', '403': '435', '339': '371', '419': '435',
                                         '323': '371', '387': '435', '371': '371', '435': '435',
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
        # convert flag to be relative to negative strand
        sam_read['FLAG'] = self.crick_mapping_conversion[sam_read['FLAG']]
        # get mapping length of the read to adjust coordinates relative to watson strand
        # convert cigar back to alpha numeric representation for mapped reads
        # change template length for paired end reads
        if sam_read['TLEN'] != '0':
            tlen = sam_read['TLEN']
            sam_read['TLEN'] = f'-{tlen}' if '-' not in tlen else tlen.replace('-', '')
