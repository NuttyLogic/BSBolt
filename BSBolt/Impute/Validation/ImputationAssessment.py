#! /usr/bin/env python3

from BSBolt.Impute.kNN_Impute import ImputeMissingValues
from BSBolt.Impute.Validation.MaskValues import MaskImputationValues


class ImputationAssessment:
    """

           :param methylation_array:
           :param masking_proportion:
           :param masking_sites:
           :param verbose:
           :param k:
           :param imputation_window_size:
           :param threads:
           :param methylation_site_order:
           :param methylation_samples:
           :param batch_size:
           :param randomize_batch:
           """

    def __init__(self, methylation_array=None, masking_proportion=0.1, masking_sites=None, verbose=False,
                 k=5, imputation_window_size=3000000, threads=4, methylation_site_order=None, methylation_samples=None,
                 batch_size=None, randomize_batch=False):
        self.methylation_array = methylation_array
        masking_kwargs = dict(methylation_array=methylation_array, masking_proportion=masking_proportion,
                              masking_sites=masking_sites, verbose=verbose)
        self.masking = MaskImputationValues(**masking_kwargs)
        self.masking_sites = masking_sites
        self.imputation_kwargs = dict(k=k, imputation_window_size=imputation_window_size, threads=threads,
                                      verbose=verbose, meth_site_order=methylation_site_order,
                                      sample_ids=methylation_samples, batch_size=batch_size,
                                      randomize_batch=randomize_batch)
        self.stats = {}
        self.imputed_array = None

    def evaluate_imputation(self):
        if self.masking_sites:
            self.masking.mask_known_sites()
        else:
            self.masking.mask_random_sites()
        self.imputation_kwargs.update(dict(meth_matrix=self.masking.methylation_array))
        self.impute_masked_values()
        self.get_imputation_statistics()

    def impute_masked_values(self):
        masked_imputation = ImputeMissingValues(**self.imputation_kwargs)
        masked_imputation.impute_values()
        self.imputed_array = masked_imputation.meth_matrix

    def get_imputation_statistics(self):
        for site, site_value in self.masking.masking_sites.items():
            row_index, column_index = (int(x) for x in site.split('_'))
            imputed_value = self.imputed_array[row_index, column_index]
            absolute_error = abs(site_value - imputed_value)
            self.stats[site] = dict(site_value=site_value, imputed_value=imputed_value, absolute_error=absolute_error)
