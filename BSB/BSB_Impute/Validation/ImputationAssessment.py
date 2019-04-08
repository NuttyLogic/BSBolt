#! /usr/bin/env python3

import numpy as np
import pandas as pd
from BSB.BSB_Impute.Imputation.GenomeImputation import GenomeImputation
from BSB.BSB_Impute.Validation.MaskValues import MaskImputationValues


class ImputationAssessment:
    """
        Keyword Arguments
            :param methylation_dataframe:
            :param masking_proportion:
            :param masking_sites:
            :param verbose:
            :param k:
            :param genome_window_size:
            :param threads:
        Attributes
            """

    def __init__(self, methylation_dataframe=None, masking_proportion=0.1, masking_sites=None, verbose=False,
                 k=5, imputation_window_size=6000000, threads=4):
        self.methylation_dataframe = methylation_dataframe
        masking_kwargs = dict(methylation_dataframe=methylation_dataframe, masking_proportion=masking_proportion,
                              masking_sites=masking_sites, verbose=verbose)
        self.masking = MaskImputationValues(**masking_kwargs)
        self.masking_sites = masking_sites
        self.imputation_kwargs = dict(k=k, imputation_window_size=imputation_window_size,
                                      threads=threads, verbose=verbose)
        self.stats = {}
        self.test_dataframe = None

    def evaluate_imputation(self):
        if self.masking_sites:
            self.masking.mask_known_sites()
        else:
            self.masking.mask_random_sites()
        evaluation_df = self.masking.methylation_dataframe
        self.imputation_kwargs.update(dict(row_labels=list(evaluation_df.index),
                                      sample_labels=list(evaluation_df),
                                      genomic_array=evaluation_df.values))
        self.impute_masked_values()
        self.get_imputation_statistics()

    def impute_masked_values(self):
        masked_imputation = GenomeImputation(**self.imputation_kwargs)
        masked_imputation.impute_windows()
        self.test_dataframe = pd.DataFrame(data=masked_imputation.genomic_array,
                                           index=list(self.masking.methylation_dataframe.index),
                                           columns=list(self.masking.methylation_dataframe))

    def get_imputation_statistics(self):
        for site, site_value in self.masking.masking_sites.items():
            row_index, column_index = (int(x) for x in site.split('_'))
            imputed_value = self.test_dataframe.iat[row_index, column_index]
            absolute_error = np.abs(site_value - imputed_value)
            self.stats[site] = dict(site_value=site_value, imputed_value=imputed_value, absolute_error=absolute_error)
