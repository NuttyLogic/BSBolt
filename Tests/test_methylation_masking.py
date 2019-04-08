import os
import unittest
import pandas as pd
from BSB.BSB_Impute.Validation.MaskValues import MaskImputationValues


test_directory = os.path.dirname(os.path.realpath(__file__))
test_methylation_data = f'{test_directory}/TestData/ch21_meth_data.tsv'
test_dataframe = pd.read_csv(test_methylation_data, sep='\t', index_col=0)

# test standard masking
masked_dataframe = MaskImputationValues(methylation_dataframe=test_dataframe, masking_proportion=0.05, verbose=True)
masked_dataframe.mask_random_sites()

# test masking known sites
new_masking = MaskImputationValues(methylation_dataframe=test_dataframe, masking_proportion=0.05, verbose=True,
                                   masking_sites=masked_dataframe.masking_sites)
new_masking.mask_known_sites()

# test masking proportions
masking_proportions = [.2] + [0 for _ in range(len(list(test_dataframe)) - 1)]
masking_proportion_test = MaskImputationValues(methylation_dataframe=test_dataframe,
                                               masking_proportion=masking_proportions,
                                               verbose=True)
masking_proportion_test.mask_random_sites()


class TestSiteMasking(unittest.TestCase):

    def setUp(self):
        pass

    def test_sites_masked(self):
        """ Test values are being masked by comparing dataframe row length before and after dropping rows with
        null values
        """
        masked_row_count = masked_dataframe.methylation_dataframe.dropna(axis=0).shape[0]
        input_row_count = test_dataframe.dropna(axis=0).shape[0]
        self.assertLess(masked_row_count, input_row_count)

    def test_same_sites_masked(self):
        """If given a list of masking sites test the same sites are masked when rerun"""
        for key in new_masking.masking_sites:
            self.assertIn(key, masked_dataframe.masking_sites)

    def test_known_value_save(self):
        """Test saved known values correspond to the data at the original index"""
        for key, value in masked_dataframe.masking_sites.items():
            row_index, column_index = (int(x) for x in key.split('_'))
            known_value = test_dataframe.iat[row_index, column_index]
            self.assertEqual(value, known_value)

    def test_masking_proportions(self):
        """Test sites masked for first site only if masking proportion provided for only that sample"""
        for key, value in masking_proportion_test.masking_sites.items():
            row_index, column_index = (int(x) for x in key.split('_'))
            self.assertEqual(column_index, 0)


if __name__ == '__main__':
    unittest.main()
