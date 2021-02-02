import os
import subprocess
import unittest

from tests.TestHelpers import test_directory


test_methylation_data = f'{test_directory}/TestData/kNN_test_matrix.txt'

imputation_command = ['python3', '-m', 'bsbolt', 'Impute', '-M', test_methylation_data, '-t', '10', '-O',
                      f'{test_directory}/TestData/test_imputed_matrix.txt', '-verbose', '-B', '3', '-R']

subprocess.run(imputation_command)


class TestCLImputation(unittest.TestCase):

    def setUp(self):
        pass

    def test_file(self):
        self.assertTrue(os.path.exists(f'{test_directory}/TestData/test_imputed_matrix.txt'))


if __name__ == '__main__':
    unittest.main()
