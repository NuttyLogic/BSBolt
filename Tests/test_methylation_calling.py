import subprocess
import time
import unittest

from BSB_CallMethylation.ProcessContigs import ProcessContigs
from BSB_CallMethylation.CallMethylation import CallMethylation
# import original calls

start = time.time()

ProcessContigs(input_file='TestData/chr22.sorted.bam',
               genome_database='TestData/ch22_test_db/hg38_chr22.fa_bowtie2',
               output_prefix='TestData/test_out',
               remove_ccgg=False,
               remove_sx_reads=False,
               text_output=True,
               min_read_depth=10,
               threads=1)

print(f'{time.time() - start:.2f}')

time.sleep(4)

test_call_methylation = CallMethylation(input_file='TestData/chr22.sorted.bam',
                                        genome_database='TestData/ch22_test_db/hg38_chr22.fa_bowtie2',
                                        remove_ccgg=False,
                                        remove_sx_reads=False,
                                        min_read_depth=10,
                                        contig='chr22')
test_call_methylation.call_methylation()


def get_cgmap_calls(file_path):
    cgmap_lines = []
    for line in open(file_path, 'r'):
        cgmap_lines.append(line.replace('\n', '').split('\t'))
    return cgmap_lines


validation_cgmap_lines = get_cgmap_calls('TestData/chr22_original_output.CGmap')
contig_test_cgmap_lines = get_cgmap_calls('TestData/test_out.CGmap')


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
