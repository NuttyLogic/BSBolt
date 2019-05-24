import unittest
import os
from BSB.BSB_Impute.Validation.ImputationAssessment import ImputationAssessment
from BSB.BSB_Impute.Impute_Utils.ImputationFunctions import get_bsb_matrix


test_directory = os.path.dirname(os.path.realpath(__file__))
test_methylation_data = f'{test_directory}/TestData/kNN_test_matrix.txt'
meth_matrix, meth_site_order, meth_samples = get_bsb_matrix(test_methylation_data)
imputation_test = ImputationAssessment(methylation_array=meth_matrix,
                                       masking_proportion=[.2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                       verbose=True,
                                       threads=8,
                                       methylation_site_order=meth_site_order,
                                       methylation_samples=meth_samples)

imputation_test.evaluate_imputation()


window_4_values = []
window_30_values = []

# manipulate values in imputation windows to calculate known neighbor

for key, value in imputation_test.stats.items():
    error = value['absolute_error']
    row, column = int(key.split('_')[0]), int(key.split('_')[1])
    if 1316 <= row <= 1425:
        window_4_values.append(error)
    elif 10066 <= row <= 11270:
        window_30_values.append(error)


class TestBatchImputation(unittest.TestCase):

    def setUp(self):
        pass

    def test_window_four(self):
        if window_4_values:
            for test_value in window_4_values:
                self.assertEqual(test_value, 0.284)

    def test_window_thirty(self):
        if window_30_values:
            for test_value in window_30_values:
                self.assertEqual(test_value, 0.4)


if __name__ == '__main__':
    unittest.main()
