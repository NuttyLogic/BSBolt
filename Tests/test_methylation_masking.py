import os
import unittest
import numpy as np
from BSB.BSB_Impute.Validation.MaskValues import MaskImputationValues
from BSB.BSB_Impute.Impute_Utils.ImputationFunctions import get_bsb_matrix


test_directory = os.path.dirname(os.path.realpath(__file__))
test_methylation_data = f'{test_directory}/TestData/kNN_test_matrix.txt'
test_matrix, test_sites, test_samples = get_bsb_matrix(test_methylation_data)

# test standard masking
random_masking = MaskImputationValues(methylation_array=test_matrix,
                                      masking_proportion=0.05, verbose=True)
random_masking.mask_random_sites()

# test masking known sites
new_masking = MaskImputationValues(methylation_array=test_matrix, masking_proportion=0.05, verbose=True,
                                   masking_sites=random_masking.masking_sites)
new_masking.mask_known_sites()

# test masking proportions
masking_proportions = [.2] + [0 for _ in range(9)]
masking_proportion_test = MaskImputationValues(methylation_array=test_matrix,
                                               masking_proportion=masking_proportions,
                                               verbose=True)
masking_proportion_test.mask_random_sites()


class TestSiteMasking(unittest.TestCase):

    def setUp(self):
        pass

    def test_sites_masked(self):
        """ Test values are being masked by comparing row length before and after dropping rows with
        null values
        """
        masked_row_count = random_masking.methylation_array[
            ~np.isnan(random_masking.methylation_array).any(axis=1)].shape[0]
        input_row_count = test_matrix[~np.isnan(test_matrix).any(axis=1)].shape[0]
        self.assertLess(masked_row_count, input_row_count)

    def test_same_sites_masked(self):
        """If given a list of masking sites test the same sites are masked when rerun"""
        for key in new_masking.masking_sites:
            self.assertIn(key, random_masking.masking_sites)

    def test_known_value_save(self):
        """Test saved known values correspond to the data at the original index"""
        for key, value in random_masking.masking_sites.items():
            row_index, column_index = (int(x) for x in key.split('_'))
            known_value = test_matrix[row_index, column_index]
            self.assertEqual(value, known_value)

    def test_masking_proportions(self):
        """Test sites masked for first site only if masking proportion provided for only that sample"""
        for key, value in masking_proportion_test.masking_sites.items():
            row_index, column_index = (int(x) for x in key.split('_'))
            self.assertEqual(column_index, 0)


if __name__ == '__main__':
    unittest.main()
