#! /usr/bin/env python3

import copy
import pickle
import unittest
from BSBolt.Align.ProcessSamReads import ProcessSamAlignment
from tests.TestHelpers import bsb_directory


# tests will only work if an index has already been generated
with open(f'{bsb_directory}tests/TestData/BSB_Test_DB/genome_index.pkl', 'rb') as lo:
    contig_lens = pickle.load(lo)


sam_read_cats = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN',
                 'SEQ', 'QUAL', 'SAM_TAGS', 'read_group', 'conversion_bases', 'mapping_reference']

read_template = {cat: None for cat in sam_read_cats}
read_template['QNAME'] = 'test_read'
read_template['SEQ'] = 'CCCCCCCTTTTTTTAAAAAAA'
read_template['QUAL'] = 'abcdefghijklmnopqrstu'
read_template['CIGAR'] = '21M'
read_template['POS'] = '10'
read_template['PNEXT'] = '20'
read_template['RNAME'] = 'chr10'
read_template['RNEXT'] = 'chr10'
read_template['TLEN'] = '10'


def get_test_read(read_flag=None, read_mismatch='XM:i:0', read_group='1',
                  read_conversion='C2T', mapping_ref='W_C2T'):
    read = copy.deepcopy(read_template)
    read['FLAG'] = read_flag
    read['SAM_TAGS'] = [read_mismatch]
    read['read_group'] = read_group
    read['conversion_bases'] = read_conversion
    read['mapping_reference'] = mapping_ref
    return read


sam_alignment_processor = ProcessSamAlignment()


class TestSamReadProcessing(unittest.TestCase):

    def setUp(self):
        pass

    def test_reverse_complement_watson(self):
        test_read = get_test_read(read_flag='147')
        sam_alignment_processor.process_read(sam_read=test_read)
        self.assertEqual(test_read['FLAG'], '131')
        self.assertEqual(test_read['SEQ'], 'TTTTTTTAAAAAAAGGGGGGG')
        self.assertEqual(test_read['QUAL'], 'abcdefghijklmnopqrstu'[::-1])

    def test_reverse_complement_crick(self):
        test_read = get_test_read(read_flag='99', mapping_ref='C_C2T')
        sam_alignment_processor.process_read(sam_read=test_read)
        self.assertEqual(test_read['FLAG'], '115')
        self.assertEqual(test_read['SEQ'], 'TTTTTTTAAAAAAAGGGGGGG')
        self.assertEqual(test_read['QUAL'], 'abcdefghijklmnopqrstu'[::-1])

    def test_sense_flag_watson(self):
        test_read = get_test_read(read_flag='99')
        sam_alignment_processor.process_read(sam_read=test_read)
        self.assertEqual(test_read['FLAG'], '67')
        self.assertEqual(test_read['SEQ'], 'CCCCCCCTTTTTTTAAAAAAA')
        self.assertEqual(test_read['QUAL'], 'abcdefghijklmnopqrstu')

    def test_anti_sense_flag_crick(self):
        test_read = get_test_read(read_flag='147', mapping_ref='C_C2T')
        test_read['CIGAR'] = '1D20M'
        sam_alignment_processor.process_read(sam_read=test_read)
        self.assertEqual(test_read['FLAG'], '179')
        self.assertEqual(test_read['SEQ'], 'CCCCCCCTTTTTTTAAAAAAA')
        self.assertEqual(test_read['QUAL'], 'abcdefghijklmnopqrstu')
        self.assertEqual(test_read['TLEN'], '-10')

    def test_mixed_unmapped_watson(self):
        test_read = get_test_read(read_flag='101')
        sam_alignment_processor.process_read(sam_read=test_read)
        self.assertEqual(test_read['FLAG'], '69')
        self.assertEqual(test_read['SEQ'], 'CCCCCCCTTTTTTTAAAAAAA')
        self.assertEqual(test_read['QUAL'], 'abcdefghijklmnopqrstu')

    def test_mixed_unmapped_crick(self):
        test_read = get_test_read(read_flag='69', mapping_ref='C_C2T')
        test_read['CIGAR'] = '1D20M'
        sam_alignment_processor.process_read(sam_read=test_read)
        self.assertEqual(test_read['FLAG'], '101')
        self.assertEqual(test_read['SEQ'], 'CCCCCCCTTTTTTTAAAAAAA')
        self.assertEqual(test_read['QUAL'], 'abcdefghijklmnopqrstu')
        # cigar str should not be modified
        self.assertEqual(test_read['CIGAR'], '1D20M')


if __name__ == '__main__':
    unittest.main()
