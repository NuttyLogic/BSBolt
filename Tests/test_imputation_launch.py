import os
import subprocess
import unittest
import pysam

# get current directory

test_directory = os.path.dirname(os.path.realpath(__file__))
bsb_directory = '/'.join(test_directory.split('/')[:-1]) + '/'
bsbolt = f'{bsb_directory}BSBolt.py'

test_methylation_data = f'{test_directory}/TestData/kNN_test_matrix.txt'

imputation_command = ['python3', bsbolt, 'Impute', '-M', test_methylation_data, '-t', '10', '-O',
                      f'{test_directory}/TestData/test_imputed_matrix.txt', '-verbose', '-B', '3', '-R']

subprocess.run(imputation_command)


class TestCLImputation(unittest.TestCase):

    def setUp(self):
        pass

    def test_file(self):
        self.assertTrue(os.path.exists(f'{test_directory}/TestData/test_imputed_matrix.txt'))


if __name__ == '__main__':
    unittest.main()
