import pickle
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

    def __init__(self, sam_line_dict=None, contig_dict=None, bsb_database=None,
                 conversion_threshold=(0.5, 5), mismatch_threshold=4):
        assert isinstance(contig_dict, dict)
        self.contig_dict = contig_dict
        self.conversion_threshold = conversion_threshold
        self.mismatch_threshold = mismatch_threshold
        self.sam_line_dict = sam_line_dict
        self.bsb_database = bsb_database
        self.read_strand_info = (('W_C2T', '+FW', False, False), ('C_C2T', '-FW', True, True),
                                 ('W_G2A', '-RC', False, False), ('C_G2A', '+RC', True, True))
        self.good_flags = {'0', '99', '147'}
        self.output_read = self.get_combined_bisulfite_read

    @property
    def get_combined_bisulfite_read(self):
        # get number of correctly maped reads
        mapped_read_number, mapped_strands = self.check_read_mapping
        # if read uniquely mapped continue
        if mapped_read_number == 1:
            # process and return read
            bisulfite_sam_read, bisulfite_strand = self.process_reads
            return mapped_read_number, bisulfite_sam_read, bisulfite_strand
        # else return mapped read number and a null bisulfite strand
        return mapped_read_number, [self.sam_line_dict['read_sequence'],
                                    self.sam_line_dict['W_C2T']['QNAME'],
                                    self.sam_line_dict['W_C2T']['QUAL'],
                                    mapped_strands], None

    @property
    def check_read_mapping(self):
        mapped_read_number = 0
        # processing assumes standard input order
        mapped_strands = []
        for read_instruction in self.read_strand_info:
            try:
                read = self.sam_line_dict[read_instruction[0]]
            except KeyError:
                break
            else:
                # if read in good flags and mismatch below threshold increase mapped_read_number
                if read['FLAG'] in self.good_flags:
                    if self.check_mismatch(read):
                        mapped_read_number += 1
                        mapped_strands.append(read_instruction[0])
        return mapped_read_number, mapped_strands

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

    @property
    def process_reads(self):
        # process reads by read group
        for read_instruction in self.read_strand_info:
            # unpack read processing instructions
            mapping_label, strand, mapping_reverse, reverse_comp = read_instruction
            # retrieve original read sequence, should always be present even is unmapped / multi-mapped
            orignal_sequence = self.sam_line_dict['read_sequence']
            try:
                # retrieve read
                read = self.sam_line_dict[read_instruction[0]]
            except KeyError:
                break
            else:
                # check for good read mapping flag
                if read['FLAG'] in self.good_flags:
                    # get contig mapping sequence
                    contig_sequence = self.get_contig_sequence(read['RNAME'])
                    # get contig length
                    contig_len = len(contig_sequence)
                    # convert cigar to numeric representation
                    cigar_tuple = convert_alpha_numeric_cigar(read['CIGAR'])
                    mapping_length = self.get_mapping_length(cigar_tuple)
                    mapping_loc = int(read['POS'])
                    quality = read['QUAL']
                    # rule to reverse flag based on mapping strand
                    if read['FLAG'] == '147':
                        if not reverse_comp:
                            reverse_comp = True
                        else:
                            reverse_comp = False
                    if reverse_comp:
                        orignal_sequence = reverse_complement(orignal_sequence)
                        # reverse cigar tuple,
                        cigar_tuple = cigar_tuple[::-1]
                        quality = quality[::-1]
                    if mapping_reverse:
                        mapping_loc = contig_len - mapping_loc - mapping_length + 2
                    mapping_genomic_sequence = contig_sequence[mapping_loc - 2: mapping_loc + len(read['SEQ'])]
                    filtered_sequence, xs, alpha_cigar = self.process_cigar_genomic_sequence(orignal_sequence,
                                                                                             mapping_genomic_sequence,
                                                                                             cigar_tuple,
                                                                                             strand)
                    read['CIGAR'] = alpha_cigar
                    read['SEQ'] = filtered_sequence
                    read['CIGAR'] = alpha_cigar
                    read['SAM_TAGS'].append(f'XS:i:{xs}')
                    read['SAM_TAGS'].append(f'XO:Z:{strand}')
                    read['POS'] = str(mapping_loc)
                    read['QUAL'] = quality
                    return dict(read), mapping_label

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

    def process_cigar_genomic_sequence(self, read_sequence, mapping_genomic_sequence, cigar, strand):
        """Take read sequence, genomic sequence and return aligned sequence according to cigar string
        Arguments:
            read_sequence (str): SAM read sequence
            mapping_genomic_sequence (str): reference sequence over mapping area
            cigar (list): list of tuples describing cigar string, can't be alpha numeric
            strand (str): reference strand for bisulfite conversion check
            """
        matched_read_sequence = []
        filtered_sequence = []
        # bisulfite conversion flag, passing by default
        xs = 0
        end_position = 0
        alpha_cigar = ''
        for cigar_type in cigar:
            # handle matching bases, extend end_position and matched_read_sequence equally
            if cigar_type[0] == 0:
                matched_read_sequence.append(read_sequence[end_position:cigar_type[1] + end_position])
                filtered_sequence.append(read_sequence[end_position:cigar_type[1] + end_position])
                end_position += cigar_type[1]
                alpha_cigar += f'{cigar_type[1]}M'
            # if insertion extend filtered sequence, and advance end position skipping insertion for matched sequence
            if cigar_type[0] == 1:
                filtered_sequence.append(read_sequence[end_position:cigar_type[1] + end_position])
                end_position += cigar_type[1]
                alpha_cigar += f'{cigar_type[1]}I'
            # if deletion append placeholder character to sequence
            if cigar_type[0] == 2:
                matched_read_sequence.append('-' * cigar_type[1])
                alpha_cigar += f'{cigar_type[1]}D'
            # if read is soft clipped extend end position, only has effect is read is reversed
            if cigar_type[0] == 4:
                end_position += cigar_type[1]
        if self.check_bisulfite_conversion(''.join(matched_read_sequence), mapping_genomic_sequence, strand):
            xs = 1
        return ''.join(filtered_sequence), xs, alpha_cigar

    def check_bisulfite_conversion(self, matched_read_sequence, mapping_genomic_sequence, strand):
        """
        Takes matched read sequence and check for bisulfite conversion of non-CpG sites in relation to the reference.
        The conversion check is strand specific, ie Watson and Crick strands will have different conversion patterns.
        Arguments
            matched_read_sequence (str): str of uppercase nucleotides with - representing deleted bases
            mapping_genomic_sequence (str): str of uppercase reference nucleotides
            strand (str): mapping strand to set conversion pattern
        Returns
            Type bool: True if nonconversion_proportion greater than threshold and the number of unconverted sites
                        great than specified threshold, else return False
        """
        # rules for + strand mapping
        ch_sequences = {'CA', 'CT', 'CC'}
        nucleotide_pattern = ('C', 'T')
        # mapping_genomic_sequence is extended by 1 bp on each end so 2bp context can be considered regardless of strand
        reference_sequence = mapping_genomic_sequence[1:]
        if strand in {'-FW', '-RC'}:
            # rules for - strand mapping
            nucleotide_pattern = ('G', 'A')
            ch_sequences = {'TG', 'AG', 'GG'}
            reference_sequence = mapping_genomic_sequence[:-1]
        converted_sites = 0.0
        unconverted_sites = 0.0
        for count in range(len(matched_read_sequence)):
            reference_nuc = matched_read_sequence[count]
            # get nucleotide context
            genomic_2mer = reference_sequence[count: count + 2]
            if genomic_2mer in ch_sequences:
                # if nucleotide is reference base, base unconverted
                if reference_nuc == nucleotide_pattern[0]:
                    unconverted_sites += 1.0
                # converted base if C > T, G > A
                elif reference_nuc == nucleotide_pattern[1]:
                    converted_sites += 1.0
        if unconverted_sites:
            nonconversion_proportion: float = unconverted_sites / (converted_sites + unconverted_sites)
        else:
            # necessary to prevent ZeroDivisionError
            nonconversion_proportion: float = 0.0
        if nonconversion_proportion >= self.conversion_threshold[0] and unconverted_sites >= self.conversion_threshold[1]:
            return True
        return False

    def get_contig_sequence(self, contig_id):
        """ Helper function to deserialize reference sequences as needed or return sequence in memory.
        Args;
            contig_id (str): contig label
        Returns:
             contig_sequence (str): str with reference sequence
        """
        try:
            contig_sequence = self.contig_dict[contig_id]
        except KeyError:
            with open(f'{self.bsb_database}{contig_id}.pkl', 'rb') as contig:
                contig_sequence = pickle.load(contig)
                self.contig_dict[contig_id] = contig_sequence
                return contig_sequence
        else:
            return contig_sequence
