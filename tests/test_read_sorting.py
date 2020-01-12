#! /usr/bin/env python3

import copy
import unittest
from BSBolt.Align.SortReads import SortReads

read_sorter = SortReads(mismatch_threshold=4, allow_discordant=False, contig_lens={'chr1': 1000})
discordant_read_sorter = SortReads(mismatch_threshold=4, allow_discordant=True, contig_lens={'chr1': 1000})
discordant_medium_mismatch = SortReads(mismatch_threshold=5, allow_discordant=True, contig_lens={'chr1': 1000})
discordant_high_mismatch = SortReads(mismatch_threshold=10, allow_discordant=True, contig_lens={'chr1': 1000})

# get current directory

sam_read_cats = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN',
                 'SEQ', 'QUAL', 'SAM_TAGS', 'read_group', 'conversion_bases', 'mapping_reference']

read_template = {cat: None for cat in sam_read_cats}
read_template['QNAME'] = 'test_read'


def get_test_read(read_flag=None, read_mismatch='XM:i:0', read_group='1',
                  read_conversion='C2T', mapping_ref='W_C2T'):
    read = copy.deepcopy(read_template)
    read['FLAG'] = read_flag
    read['SAM_TAGS'] = [read_mismatch]
    read['read_group'] = read_group
    read['conversion_bases'] = read_conversion
    read['mapping_reference'] = mapping_ref
    read['CIGAR'] = '50M'
    read['POS'] = '100'
    read['RNAME'] = 'chr1'
    return read

# return unmapped pared end reads if discordant not allowed, mixed if allowed


paired_reads_1 = []

pe_flags = ['99', '147', '83', '163']
pe_mismatches = ['XM:i:0', 'XM:i:5', 'XM:i:10', 'XM:i:10']
pe_conversions = ['C2T', 'G2A', 'G2A', 'C2T']
pe_read_groups = ['1', '1', '2', '2']
pe_mapping_references = ['W_C2T', 'W_C2T', 'C_G2A', 'C_G2A']
pe_mapping_settings = [pe_flags, pe_mismatches, pe_conversions, pe_read_groups, pe_mapping_references]

for flag, mismatch, conversion, group, ref in zip(*pe_mapping_settings):
    paired_reads_1.append(get_test_read(read_flag=flag, read_mismatch=mismatch,
                                        read_conversion=conversion, read_group=group,
                                        mapping_ref=ref))

discordant_paired_reads = []

pe_dis_flags = ['99', '147', '355', '403']
pe_dis_mismatches = ['XM:i:0', 'XM:i:5', 'XM:i:0', 'XM:i:5']
pe_dis_conversions = ['C2T', 'G2A', 'G2A', 'C2T']
pe_dis_read_groups = ['1', '1', '1', '1']
pe_dis_mapping_references = ['W_C2T', 'W_C2T', 'W_C2T', 'W_C2T']
pe_dis_mapping_settings = [pe_dis_flags, pe_dis_mismatches, pe_dis_conversions,
                           pe_dis_read_groups, pe_dis_mapping_references]

for flag, mismatch, conversion, group, ref in zip(*pe_dis_mapping_settings):
    discordant_paired_reads.append(get_test_read(read_flag=flag, read_mismatch=mismatch,
                                                 read_conversion=conversion, read_group=group,
                                                 mapping_ref=ref))


single_end_mapped_reads = []

single_flags = ['0', '256', '0']
single_mismatch = ['XM:i:0', 'XM:i:0', 'XM:i:5']
single_conversion = ['C2T', 'C2T', 'G2A']
single_mapping_ref = ['C_C2T', 'C_C2T', 'W_G2A']
single_read_group = ['1', '1', '2']
single_mapping_stats = [single_flags, single_mismatch, single_conversion,
                        single_mapping_ref, single_read_group]

for flag, mismatch, conversion, ref, group in zip(*single_mapping_stats):
    single_end_mapped_reads.append(get_test_read(read_flag=flag,
                                                 read_mismatch=mismatch,
                                                 read_conversion=conversion,
                                                 mapping_ref=ref,
                                                 read_group=group))


class TestReadSorting(unittest.TestCase):

    def setUp(self):
        pass

    def test_non_discordant_no_mapping_pair(self):
        # expect unmapped read pair with flags set to non-mapping
        non_discordant_pair = read_sorter.get_reads(read_group=copy.deepcopy(paired_reads_1))
        # no mapping reads
        self.assertFalse(non_discordant_pair[0])
        # assert no mapping_references
        self.assertIsNone(non_discordant_pair[1])
        # assert flags are properly set to non-mapping_flags
        self.assertEqual(non_discordant_pair[2][0]['FLAG'], '77')
        self.assertEqual(non_discordant_pair[2][1]['FLAG'], '141')

    def test_mixed_mapping_pair(self):
        # read return as mixed pair with first read mapped to sense strand
        mixed_pair = discordant_read_sorter.get_reads(read_group=copy.deepcopy(paired_reads_1))
        # assert mixed mapping
        self.assertEqual(mixed_pair[0], 'Mixed1')
        # assert no mapping ref returned
        self.assertIsNone(mixed_pair[1])
        # assert flags are properly set
        self.assertEqual(mixed_pair[2][0]['FLAG'], '73')
        self.assertEqual(mixed_pair[2][1]['FLAG'], '133')

    def test_proper_mapping_pair(self):
        # read proper read pairing
        proper_pair = discordant_medium_mismatch.get_reads(read_group=copy.deepcopy(paired_reads_1))
        # assert mixed mapping
        self.assertEqual(proper_pair[0], 'Paired')
        # assert no mapping ref returned
        self.assertIn('W_C2T', proper_pair[1])
        # assert flags are properly set
        self.assertEqual(proper_pair[2][0]['FLAG'], '99')
        self.assertEqual(proper_pair[2][1]['FLAG'], '147')

    def test_unmapped_multireference_pair(self):
        # read pairs map to multiple references
        unmapped_multireference = discordant_high_mismatch.get_reads(read_group=copy.deepcopy(paired_reads_1))
        # assert mapping is false
        self.assertFalse(unmapped_multireference[0])
        # mapping ref should contain both references
        self.assertIn('W_C2T', unmapped_multireference[1])
        self.assertIn('C_G2A', unmapped_multireference[1])
        # assert flags are unmapped
        self.assertEqual(unmapped_multireference[2][0]['FLAG'], '77')
        self.assertEqual(unmapped_multireference[2][1]['FLAG'], '141')

    def test_discordant_read_pairing(self):
        discordant_read_pairs = discordant_read_sorter.get_reads(read_group=copy.deepcopy(discordant_paired_reads))
        # assert reads are discordant
        self.assertEqual(discordant_read_pairs[0], 'Discordant')
        # assert not mapping ref
        self.assertIsNone(discordant_read_pairs[1])
        # assert expected flags appear
        self.assertEqual(discordant_read_pairs[2][0]['FLAG'], '65')
        self.assertEqual(discordant_read_pairs[2][1]['FLAG'], '129')

    def test_discordant_proper_pairing(self):
        discordant_proper = discordant_medium_mismatch.get_reads(read_group=copy.deepcopy(discordant_paired_reads))
        # assert paired mapping
        self.assertEqual(discordant_proper[0], 'Paired')
        # assert mapping ref
        self.assertIn('W_C2T', discordant_proper[1])
        # assert flags
        self.assertEqual(discordant_proper[2][0]['FLAG'], '99')
        self.assertEqual(discordant_proper[2][1]['FLAG'], '147')
        self.assertEqual(discordant_proper[2][2]['FLAG'], '355')
        self.assertEqual(discordant_proper[2][3]['FLAG'], '403')

    def test_single_mapping(self):
        single_mapping = read_sorter.get_reads(read_group=copy.deepcopy(single_end_mapped_reads))
        # assert paired mapping
        self.assertEqual(single_mapping[0], 'Single')
        # assert mapping ref
        self.assertIn('C_C2T', single_mapping[1])
        # assert flags
        self.assertEqual(single_mapping[2][0]['FLAG'], '0')
        self.assertEqual(single_mapping[2][1]['FLAG'], '256')

    def test_single_multireference(self):
        single_multi = discordant_medium_mismatch.get_reads(read_group=copy.deepcopy(single_end_mapped_reads))
        # assert no mapping
        self.assertFalse(single_multi[0])
        # assert two mapping ref
        self.assertIn('C_C2T', single_multi[1])
        self.assertIn('W_G2A', single_multi[1])
        # assert unmapped flag
        self.assertEqual(single_multi[2][0]['FLAG'], '4')


if __name__ == '__main__':
    unittest.main()
