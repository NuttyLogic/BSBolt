from BSB.BSB_Align.AlignmentHelpers import convert_alpha_numeric_cigar
from BSB.BSB_Utils.UtilityFunctions import reverse_complement


class ProcessSamAlignment:
    """ SAM alignment processing class. Takes combined fastq and sam reads. Outputs process SAM reads, if read is
    uniquely aligned the correctly mapped read will be sorted. If read is multimapped or unmapped the original fastq
    sequence will be output.
    Keyword Arguments:
            sam_line_dict (dict): dict of import sam and fastq reads
            contig_dict (dict): dict of contig sequences, version of reference genome loaded in memory
            bsb_database (str): path to directory where bsb database generated
            conversion_threshold (tuple): conversion proprotion, conversion threshold to consider an incompletely
                                        converted read
            mismatch_threshold (int): number of mismatches allowed for a read to be considered
    Attributes:
        self.sam_line_dict (dict): dict of import sam and fastq reads
        self.contig_dict (dict): dict of contig sequences, version of reference genome loaded in memory
        self.bsb_database (str): path to directory where bsb database generated
        self.conversion_threshold (tuple): conversion proprotion, conversion threshold to consider an incompletely
                                    converted read
       self.mismatch_threshold (int): number of mismatches allowed for a read to be considered
       self.read_strand_info (tuple): Stored information about strand specific read processing
       self.good_flags (set): sam flags that indicate a properly aligned read or read pair
       self.output_read (object): processed sam read handled by BSB_Align

    """

    def __init__(self, sam_reads=None, contig_lens=None, mismatch_threshold=4):
        self.sam_reads = sam_reads
        self.mismatch_threshold = mismatch_threshold
        self.contig_lens = contig_lens
        self.mapping_flags = {'0', '16', '83', '99', '147', '163'}
        self.crick_mapping_conversion = {'0': '16', '16': '0',
                                         '83': '99', '163': '147',
                                         '99': '83', '147': '163'}
        self.output_reads = self.get_combined_bisulfite_read

    @property
    def get_combined_bisulfite_read(self):
        # get number of correctly maped reads
        mapped_read_number, mapped_reads = self.check_read_mapping
        if mapped_read_number == 1:
            processed_reads = []
            for read_grouping in mapped_reads:
                reads = [self.process_read(read, mapped_read_number) for read in read_grouping]
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
                    if not self.check_mismatch(read):
                        mismatch_check = True
                if not mismatch_check:
                    mapped_strands.add(read_group[0]['mapping_reference'])
                    mapping_reads.append(read_group)
        if len(mapped_strands) == 0:
            return len(mapped_strands), [self.sam_reads[0]]
        elif len(mapped_strands) == 1:
            return len(mapped_strands), mapping_reads
        return len(mapped_strands), self.set_unmapped_read(self.sam_reads[0], mapped_strands)

    @staticmethod
    def set_unmapped_read(read_group, mapping_strands):
        unmapped_flags = tuple('4') if len(read_group) == 1 else ('77', '141')
        mapping_stand_id = ','.join(mapping_strands)
        unmapped_reads = []
        for read, flag in zip(read_group, unmapped_flags):
            unmapped_read = dict(read)
            unmapped_read['FLAG'] = flag
            unmapped_read['RNAME'] = '*'
            unmapped_read['POS'] = 0
            unmapped_read['CIGAR'] = 0
            unmapped_read['PNEXT'] = 0
            unmapped_read['TLEN'] = 0
            unmapped_reads['SAM_TAGS'] = ['YT:Z:UP', f'XO:Z:{mapping_stand_id}']
            unmapped_reads.append(unmapped_read)
        return tuple(unmapped_reads)

    def check_mismatch(self, read):
        """ Arguments:
                read (dict): individual sam read
        """
        # retrieve sam tags
        SAM_TAGS = read['SAM_TAGS']
        for tag in SAM_TAGS:
            # if sam tag XM, mismatch sam tag
            if 'XM:i' in tag:
                # parse sam tag
                mismatch_number = int(tag.split(':')[-1])
                if mismatch_number <= self.mismatch_threshold:
                    return True
                # set bad flag and return false
                read['FLAG'] = '4'
                return False
        # if mismatch sam flag missing return false by default
        return False

    def process_read(self, read, paired_end):
        mapped_read = dict(read)
        contig_len = self.contig_lens[mapped_read['RNAME']]
        sequence = mapped_read['SEQ']
        cigar_tuple = convert_alpha_numeric_cigar(mapped_read['CIGAR'])
        mapping_length = self.get_mapping_length(cigar_tuple)
        mapping_loc = int(mapped_read['POS'])
        quality = mapped_read['QUAL']
        # rule to reverse flag based on mapping strand
        if 'C_' in mapped_read['mapping_reference']:
            sequence = reverse_complement(sequence)
            # reverse cigar tuple,
            cigar_tuple = cigar_tuple[::-1]
            quality = quality[::-1]
            mapping_loc = contig_len - mapping_loc - mapping_length + 2
            if paired_end > 1:
                tlen = mapped_read['TLEN']
                mapped_read['TLEN'] = tlen.replace('-', '') if '-' in tlen else f'-{tlen}'
        filtered_sequence, alpha_cigar, filtered_qual = self.process_cigar_genomic_sequence(sequence,
                                                                                            cigar_tuple,
                                                                                            quality)
        read['CIGAR'] = alpha_cigar
        read['SEQ'] = filtered_sequence
        read['CIGAR'] = alpha_cigar
        read['SAM_TAGS'].append(f'XO:Z:{mapped_read["mapping_reference"]}')
        read['POS'] = str(mapping_loc)
        read['QUAL'] = filtered_qual
        return mapped_read

    @staticmethod
    def get_mapping_length(cigar_tuple):
        """Return effective mapping length, ie. the representation of read to properly place reverse complemented reads
         """
        mapping_length = 0
        # if sam flag is match or deletion add to mapping length
        counting_cigar_labels = {0, 2}
        for cigar_label in cigar_tuple:
            if cigar_label[0] in counting_cigar_labels:
                mapping_length += cigar_label[1]
        return mapping_length

    def process_cigar_genomic_sequence(self, read_sequence, cigar, qual):
        """Take read sequence, genomic sequence and return aligned sequence according to cigar string
        Arguments:
            read_sequence (str): SAM read sequence
            cigar (list): list of tuples describing cigar string, can't be alpha numeric
            strand (str): reference strand for bisulfite conversion check
            """
        matched_read_sequence = []
        filtered_sequence = []
        filtered_qual = []
        end_position = 0
        alpha_cigar = ''
        for cigar_type in cigar:
            # handle matching bases, extend end_position and matched_read_sequence equally
            if cigar_type[0] == 0:
                matched_read_sequence.append(read_sequence[end_position:cigar_type[1] + end_position])
                filtered_sequence.append(read_sequence[end_position:cigar_type[1] + end_position])
                filtered_qual.append(qual[end_position:cigar_type[1] + end_position])
                end_position += cigar_type[1]
                alpha_cigar += f'{cigar_type[1]}M'
            # if insertion extend filtered sequence, and advance end position skipping insertion for matched sequence
            if cigar_type[0] == 1:
                filtered_sequence.append(read_sequence[end_position:cigar_type[1] + end_position])
                filtered_qual.append(qual[end_position:cigar_type[1] + end_position])
                end_position += cigar_type[1]
                alpha_cigar += f'{cigar_type[1]}I'
            # if deletion append placeholder character to sequence
            if cigar_type[0] == 2:
                matched_read_sequence.append('-' * cigar_type[1])
                alpha_cigar += f'{cigar_type[1]}D'
            # if read is soft clipped extend end position, only has effect is read is reversed
            if cigar_type[0] == 4:
                end_position += cigar_type[1]
        filtered_sequence, filtered_qual = ''.join(filtered_sequence), ''.join(filtered_qual)
        assert len(filtered_sequence) == len(filtered_qual)
        return filtered_sequence, alpha_cigar, filtered_qual
