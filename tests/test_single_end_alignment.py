import os
import subprocess
import unittest

# get current directory

test_directory = os.path.dirname(os.path.realpath(__file__))
bsb_directory = '/'.join(test_directory.split('/')[:-1]) + '/'
bsbolt = f'{bsb_directory}BSBolt.py'
# generate simulated reads


bsb_align_commands = ['python3', bsbolt, 'Align',
                      '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB', '-F1',
                      f'{bsb_directory}Tests/TestSimulations/BSB_pe_meth_1.fastq', '-O',
                      f'{bsb_directory}Tests/BSB_pe_test', '-S', '-OU', '-BT2-k', '10', '-BT2-p', '10']
subprocess.run(bsb_align_commands)


class TestSingleEndAlignment(unittest.TestCase):

    def setup(self):
        pass

    def test_file_exists(self):
        self.assertTrue(os.path.isfile(f'{bsb_directory}Tests/BSB_pe_test.sorted.bam'))


if __name__ == '__main__':
    unittest.main()
