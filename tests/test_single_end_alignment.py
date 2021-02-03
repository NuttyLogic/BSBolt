import os
import subprocess
import unittest
from tests.TestHelpers import bsb_directory


bsb_align_commands = ['python3', '-m', 'bsbolt', 'Align',
                      '-DB', f'{bsb_directory}tests/TestData/BSB_Test_DB', '-F1',
                      f'{bsb_directory}tests/TestSimulations/BSB_pe_1.fq', '-O',
                      f'{bsb_directory}tests/BSB_se_test', '-t', '10', '-UN']
subprocess.run(bsb_align_commands)


class TestSingleEndAlignment(unittest.TestCase):

    def setup(self):
        pass

    def test_file_exists(self):
        self.assertTrue(os.path.exists(f'{bsb_directory}tests/BSB_se_test.bam'))


if __name__ == '__main__':
    unittest.main()
