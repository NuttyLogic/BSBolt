import unittest
import os
import pandas as pd
from BSB.BSB_Impute.kNN_Impute import ImputeMissingValues
from BSB.BSB_Impute.Validation.MaskValues import MaskImputationValues

test_directory = os.path.dirname(os.path.realpath(__file__))
test_methylation_data = f'{test_directory}/TestData/ch21_meth_data.tsv'
test_df = pd.read_csv(test_methylation_data, sep='\t', index_col=0)
masking = MaskImputationValues(methylation_dataframe=test_df, masking_proportion=.1)
masking.mask_random_sites()
masked_df = masking.methylation_dataframe
print(test_df.values)
print(masked_df.dropna(axis=0).shape)
print(masked_df.values)
test_imputation = ImputeMissingValues(batch_size=50)
test_imputation.meth_matrix = masked_df.values
test_imputation.meth_site_order = list(masked_df.index)
test_imputation.sample_ids = ['blah', list(masked_df)]
test_imputation.impute_values()
print(masked_df.dropna(axis=0).shape)
print(masked_df.values)
print(test_imputation.meth_matrix)





class TestBatchImputation(unittest.TestCase):

    def setUp(self):
        pass

