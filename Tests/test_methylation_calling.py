import unittest
import pysam
from BSB.BSB_CallMethylation.CallMethylation import CallMethylation


class TestMethylationCalling(unittest.TestCase):

    def setUp(self):
        pass

    def test_context_tables(self):
        test_context = test_call_methylation.context_tables['context_table']['CAA']
        self.assertEqual(test_context, 'CHH')
        test_subcontext = test_call_methylation.context_tables['sub_context_table']['CAC']
        self.assertEqual(test_subcontext, 'CA')

    def test_cgmap_line(self):
        for validation_line, contig_test_line in zip(validation_cgmap_lines, contig_test_cgmap_lines):
            self.assertEqual(validation_line[2], contig_test_line[2])
            self.assertAlmostEqual(float(validation_line[5]), float(contig_test_line[5]), places=1)


if __name__ == '__main__':
    unittest.main()
