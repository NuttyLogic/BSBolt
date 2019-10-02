from BSB.BSB_Align.SortDiscordantReads import SortDiscordantReads


class ReadSortingError(Exception):
    """Error in methylation calling process"""
    pass


class SortReads(object):

    def __init__(self, mismatch_threshold=4, allow_discordant=True):
        self.mismatch_threshold = mismatch_threshold
        self.allow_discordant = allow_discordant
        self.discordant_sorter = SortDiscordantReads()

    def get_reads(self, read_group):
        mapping_references, read_copies = self.check_mapping_status(read_group=read_group)
        # return single ended mapped or unmapped reads
        if len(read_copies) == 1:
            if not mapping_references:
                return False, None, self.set_unmapped_read(read_copies[0],
                                                           mapping_strands=list(mapping_references))
            return self.handle_single_end(mapping_references, read_copies)
        # if no mapped paired end reads return unmapped reads
        if not mapping_references:
            return False, mapping_references, self.set_unmapped_read(read_copies,
                                                                     mapping_strands=list(mapping_references))
        # look for properly paired reads
        mapping_references, paired_reads = self.check_proper_read_pairing(read_group)
        if len(mapping_references) > 1:
            return False, mapping_references, self.set_unmapped_read(read_copies,
                                                                     mapping_strands=list(mapping_references))
        elif len(mapping_references) == 1:
            return 'Paired', mapping_references, paired_reads
        # check for mixed and discordant reads, start by only considering valid reads
        read_group = [read for read in read_group if read['mapping_status']]
        if not read_group or not self.allow_discordant:
            return False, None, self.set_unmapped_read(read_copies, mapping_strands=[])
        discord_status, discord_mapping, discordant_reads = self.check_discordant(read_group, read_copies)
        if not discordant_reads:
            return False, discord_mapping, self.set_unmapped_read(read_copies, mapping_strands=discord_mapping)
        return discord_status, discord_mapping, discordant_reads

    def handle_single_end(self, mapping_references, read_copies, read_group):
        if len(mapping_references) != 1:
            return False, mapping_references, self.set_unmapped_read(read_copies[0],
                                                                     mapping_strands=list(mapping_references))
        else:
            return 'Single', mapping_references, [read for read in read_group if read['mapping_status']]

    def check_mapping_status(self, read_group):
        read_1, read_2 = read_group[0], None
        mate_conversion = 'G2A' if read_1['conversion_bases'] == 'C2T' else 'C2T'
        mapping_references = set()
        for read in read_group:
            if read['read_group'] == '1' and read['conversion_bases'] == mate_conversion:
                read_2 = read
            if self.check_mismatch(read):
                read['mapping_status'] = False
            else:
                read['mapping_status'] = True
                mapping_references.add(read['mapping_reference'])
        read_copies = [dict(read_1)]
        if read_2:
            read_copies.append(dict(read_2))
        return mapping_references, read_copies

    @staticmethod
    def check_proper_read_pairing(read_group):
        """
        Return sorted proper reads pair if present.
        Prefer valid read pairs to discordant / mixed read mappings
        """
        properly_paired_flags = {'83', '99', '147', '163', '256', '272', '339', '355', '403', '419'}
        paired_reads = []
        for read in read_group:
            if read['FLAG'] in properly_paired_flags:
                paired_reads.append(read)
        mapping_references = set()
        sorted_reads = []
        for read_1, read_2 in zip(paired_reads[::2], paired_reads[1::2]):
            if read_1['mapping_status'] and read_2['mapping_status']:
                assert read_1['mapping_reference'] == read_2['mapping_reference']
                mapping_references.add(read_1['mapping_reference'])
                sorted_reads.extend([read_1, read_2])
        return mapping_references, sorted_reads

    def check_discordant(self, read_group, read_copies):
        """Return discordant reads, read that don't map with correct inter-mate distance or strand orientation is
        incorrect, but both reads map uniquely. All read should come with a matched pair"""
        return self.discordant_sorter.check_discordant_read_pairing(read_group=read_group, read_copies=read_copies)

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

    @staticmethod
    def set_unmapped_read(read_group, mapping_strands):
        unmapped_read_flags = {'4', '77', '141'}
        unmapped_flags = tuple('4') if len(read_group) == 1 else ('77', '141')
        mapping_strand_id = ','.join(mapping_strands)
        unmapped_reads = []
        for read, flag in zip(read_group, unmapped_flags):
            if read['FLAG'] in unmapped_read_flags:
                unmapped_reads.append(read)
            else:
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
        return unmapped_reads
