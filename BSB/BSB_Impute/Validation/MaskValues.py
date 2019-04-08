import random
import pandas as pd
import numpy as np
from tqdm import tqdm


class MaskImputationValues:
    """
    Mask known beta values from input pd.DataFrame to evaluate imputation accuracy
        Keyword Arguments
            methylation_dataframe (pd.DataFrame): input pd.DataFrame of methylation values
            masking_proportion (float, list): proportion of sites to mask, per sample or overall
            masking_sites (dict): dict of sites that were masked hashing to known value
            verbose (bool): use tqdm progress bar or suppress
        Attributes:
            self.methylation_dataframe (pd.DataFrame): input pd.DataFrame of methylation values
            self.masking_proportion (float, list): proportion of sites to mask, per sample or overall
            self.masking_sites (dict): dict of sites that were masked hashing to known value
            self.disable_tqdm (bool): if verbose False else True
    """

    def __init__(self, methylation_dataframe=None, masking_proportion=None, masking_sites=None, verbose=False):
        assert isinstance(methylation_dataframe, pd.DataFrame)
        assert methylation_dataframe.shape[0] > methylation_dataframe.shape[1], \
            f'Number of Genomic Locations {methylation_dataframe.shape[0]} < Number of Samples ' \
            f'{methylation_dataframe.shape[1]}, Transpose Matrix if Rows Represent Samples instead and not ' \
            f'Genomic Locations'
        self.methylation_dataframe = methylation_dataframe.copy(deep=True)
        self.masking_proportion = masking_proportion
        self.methylation_samples = list(methylation_dataframe)
        self.masking_sites = masking_sites
        if self.masking_sites:
            self.masking_sites = dict(self.masking_sites)
        self.tqdm_disable = True
        if verbose:
            self.tqdm_disable = False

    def get_masking_proportion(self, sample_index):
        """ Return masking probability for site, return self.masking_probabiliyt if float given else return
        item at self.masking_proportion[sample_index]"""
        if isinstance(self.masking_proportion, float):
            return self.masking_proportion
        else:
            assert len(self.methylation_samples) == len(self.masking_proportion), 'Masking proportion mush be ' \
                                                                                  'provided for every sample'
            return self.masking_proportion[sample_index]

    @staticmethod
    def mask_value(masking_proportion):
        """Return True if random.random() less than masking proportion"""
        return random.random() < masking_proportion

    def mask_random_sites(self):
        """Mask random sites if called based on masking_proportion provided"""
        # ensure self.masking_sites is empty if called
        self.masking_sites = {}
        row_index = 0
        # iterate through rows
        for row, row_items in tqdm(self.methylation_dataframe.iterrows(), total=self.methylation_dataframe.shape[0],
                                   disable=self.tqdm_disable, desc='Masking Sites'):
            row_values = row_items.values
            for sample_index, sample_value in enumerate(row_values):
                masking_proportion = self.get_masking_proportion(sample_index)
                # if true mask site
                if self.mask_value(masking_proportion):
                    # update known value
                    self.masking_sites[f'{row_index}_{sample_index}'] = sample_value
                    row_values[sample_index] = np.nan
            row_index += 1

    def mask_known_sites(self):
        """Mask sites using provided dict of previously masked sites"""
        row_index = 0
        for row, row_items in tqdm(self.methylation_dataframe.iterrows(), total=self.methylation_dataframe.shape[0],
                                   disable=self.tqdm_disable, desc='Masking Sites'):
            row_values = row_items.values
            for sample_index, sample_value in enumerate(row_values):
                if f'{row_index}_{sample_index}' in self.masking_sites:
                    row_values[sample_index] = np.nan
            row_index += 1
