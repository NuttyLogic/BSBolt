class SortDiscordantReads:
    """
    Sorts mapped reads that pass quality control filters, but a mate pair is filtered or reads didn't map as a pair
    originally. Discordant reads, then mixed reads are returned.
    Attributes:
        self.sense_flags (set): flags of reads mapped to sense strand
        self.non_primary_flag_conversion (dict): discordant flags plus 256 to indicate non-primary alignment
    """

    def __init__(self):
        self.sense_flags = {'65', '67', '73', '97', '99', '129', '131', '137', '161', '163',
                            '321', '323', '329', '353', '355', '385', '387', '393', '417', '419'}
        self.non_primary_flag_conversion = {'65': '321', '129': '385', '97': '353', '145': '401',
                                            '81': '337', '161': '417', '113': '369', '177': '433',
                                            '73': '329', '137': '393', '89': '345', '153': '409'}

    def check_discordant_read_pairing(self, read_group, read_copies):
        """Return discordant reads, read that don't map with correct inter-mate distance or strand orientation is
        incorrect, but both reads map uniquely. All read should come with a matched pair"""
        first_reads, second_reads = [], []
        first_conversion_bases = read_copies[0]['conversion_bases']
        for read in read_group:
            if read['read_group'] == 2:
                first_conversion_bases = read_copies[1]['conversion_bases']
            if read['conversion_bases'] == first_conversion_bases:
                first_reads.append(read)
            else:
                second_reads.append(read)
        first_mapping, first_reference = self.check_discordant_pair_read_group(first_reads)
        second_mapping, second_references = self.check_discordant_pair_read_group(second_reads)
        # if both first and second mappings reads are discordant
        if not first_mapping and not second_mapping:
            combined_ref = list(set(first_reference + second_references))
            return None, combined_ref, None
        elif first_mapping and second_mapping:
            return 'Discordant', None, self.get_discordant_reads(first_reads, second_reads)
        else:
            if first_mapping:
                return 'Mixed1', None, self.get_mixed_reads(first_reads, read_copies, first=True)
            else:
                return 'Mixed2', None, self.get_mixed_reads(second_reads, read_copies, first=False)

    @staticmethod
    def check_discordant_pair_read_group(pair_read_group):
        """Reads within a group can only have one mapping reference"""
        mapping_references = set()
        for read in pair_read_group:
            mapping_references.add(read['mapping_reference'])
        if len(mapping_references) == 1:
            return True, mapping_references.pop()
        return False, list(mapping_references)

    def get_discordant_reads(self, first_reads, second_reads):
        primary_first_read, primary_second_read = first_reads[0], second_reads[0]
        primary_flags = self.get_discordant_flag(primary_first_read, primary_second_read)
        primary_first_read['FLAG'], primary_second_read['FLAG'] = primary_flags
        self.set_dicordant_read(primary_first_read, primary_second_read)
        self.set_dicordant_read(primary_second_read, primary_first_read)
        for read in first_reads[1:]:
            self.set_dicordant_read(read, primary_second_read)
            flag, _ = self.get_discordant_flag(read, primary_second_read)
            read['FLAG'] = self.non_primary_flag_conversion[flag]
        for read in second_reads[1:]:
            self.set_dicordant_read(read, primary_first_read)
            _, flag = self.get_discordant_flag(primary_first_read, read)
            read['FLAG'] = self.non_primary_flag_conversion[flag]
        return first_reads + second_reads

    @staticmethod
    def set_dicordant_read(read, reference_read):
        if read['RNAME'] == reference_read['RNAME']:
            read['RNEXT'] = '='
        else:
            read['RNEXT'] = str(reference_read['RNAME'])
        read['PNEXT'] = str(reference_read['POS'])
        read['TLEN'] = '0'

    def get_discordant_flag(self, first_read, second_read):
        first_strand = first_read['FLAG'] in self.sense_flags
        second_strand = second_read['FLAG'] in self.sense_flags
        # both reads mapped to sense strand
        if first_strand and second_strand:
            return '65', '129'
        elif first_strand and not second_strand:
            return '97', '145'
        elif not first_strand and second_strand:
            return '81', '161'
        else:
            return '113', '177'

    def get_mixed_reads(self, read_group, read_copies, first=True):
        unmapped_read = read_copies[1] if first else read_copies[0]
        primary_mapped = read_group[0]
        mapped_flag, unmapped_flag = self.get_mixed_flags(primary_mapped, first=first)
        primary_mapped['FLAG'], unmapped_read['FLAG'] = mapped_flag, unmapped_flag
        self.set_unmapped_mixed(primary_mapped, unmapped_read)
        self.set_mapped_mixed(primary_mapped)
        for read in read_group[1:]:
            read['FLAG'] = self.non_primary_flag_conversion[self.get_mixed_flags(read, first=first)[0]]
            self.set_mapped_mixed(read)
        return read_group + [unmapped_read]

    @staticmethod
    def set_mapped_mixed(mapped_read):
        mapped_read['PNEXT'] = str(mapped_read['POS'])
        mapped_read['RNEXT'] = '='
        mapped_read['TLEN'] = '0'

    @staticmethod
    def set_unmapped_mixed(mapped_read, unmapped_read):
        unmapped_read['MAPQ'] = '0'
        unmapped_read['CIGAR'] = '*'
        unmapped_read['POS'] = str(mapped_read['POS'])
        unmapped_read['PNEXT'] = str(mapped_read['POS'])
        unmapped_read['RNAME'] = str(mapped_read['RNAME'])
        unmapped_read['RNEXT'] = '='
        unmapped_read['SAM_TAGS'] = ['YT:Z:UP']
        unmapped_read['TLEN'] = '0'
        unmapped_read['mapping_reference'] = str(mapped_read['mapping_reference'])

    def get_mixed_flags(self, read, first=True):
        if read['FLAG'] in self.sense_flags:
            if first:
                return '73', '133'
            else:
                return '137', '69'
        else:
            if first:
                return '89', '165'
            else:
                return '153', '101'
