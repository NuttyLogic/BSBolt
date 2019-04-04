import unittest
import os
from BSB.BSB_Impute.kNN_Impute import ImputeMissingValues

test_directory = os.path.dirname(os.path.realpath(__file__))
test_methylation_data = f'{test_directory}/TestData/ch21_meth_data.tsv'


class TestBatchImputation(unittest.TestCase):

    def setUp(self):
        pass

    def test_sample_labels(self):
        test_batch_imputation = ImputeMissingValues(methylation_data=test_methylation_data,
                                                    output_path=f'{test_methylation_data}_test_imputation',
                                                    batch_size=10)
        self.assertEqual(len(test_batch_imputation.sample_labels), 109)

    def test_batch_size_large_remainder(self):
        test_batch_imputation = ImputeMissingValues(methylation_data=test_methylation_data,
                                                    output_path=f'{test_methylation_data}_test_imputation',
                                                    batch_size=10)
        test_batches = test_batch_imputation.get_batch_split
        self.assertEqual(test_batches[-1], 9)

    def test_batch_size_small_remainder(self):
        test_batch_imputation = ImputeMissingValues(methylation_data=test_methylation_data,
                                                    output_path=f'{test_methylation_data}_test_imputation',
                                                    batch_size=20)
        test_batches = test_batch_imputation.get_batch_split
        self.assertEqual(test_batches[0], 21)
        self.assertEqual(test_batches[-1], 25)

    def test_batch_imputation(self):
        test_batch_imputation = ImputeMissingValues(methylation_data=test_methylation_data,
                                                    output_path=f'{test_methylation_data}_test_imputation',
                                                    batch_size=20)
        test_batch_imputation.impute_values()


if __name__ == '__main__':
    unittest.main()
