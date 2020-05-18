#! /user/bin/env python3

import multiprocessing as mp
from typing import List
import numpy as np
from tqdm import tqdm
from BSBolt.Impute.Imputation.GenomeImputationWindows import GenomeImputationWindows
from BSBolt.Impute.Imputation.EuclideanDistance import get_euclidean


class ImputationError(Exception):
    """Error in imputation process"""
    pass


def imputation_process_error(error):
    """Raise if exception thrown in methylation calling process"""
    raise ImputationError(error)


class GenomeImputation:
    """kNN imputation wrapper. Handles calculating nearest neighbors and setting missing values.
       KeywordArguments:
        row_labels (list): CpG Sites structured as chrom:pos
        sample_labels (list): sample labels
        genomic_array (np.array): array of methylation values with np.NaN representing empty sites
        imputation_window_size (int): size of total distance window, values are imputed for the middle third
        k (int): number of neighbors to pull values from
        threads (int): number of threads to calculate pairwise distance
        verbose (bool): default=False, tqdm output
       Attributes
        self.row_labels (list): CpG Sites structured as chrom:pos
        self.sample_labels (list): sample labels
        self.genomic_array (np.array): array of methylation values with np.NaN representing empty sites
        self.imputation_window_size (int): size of total distance window, values are imputed for the middle third
        self.k (int): number of neighbors to pull values from
        self.threads (int): number of threads to calculate pairwise distance
        self.verbose (bool): default=False, tqdm output
        self.pool (mp.pool): multiprocessing kNN pool
        self.global_neighbors (dict): reference values to use for empty imputation windows
        self.neighbors_dicts (dict) imputation pairwise distance neighbors
        self.genomic_array_t (np.array): reference to transposed version of self.genomic_array
           """

    def __init__(self, row_labels: List[str] = None, sample_labels: List[str] = None,
                 genomic_array: np.ndarray = None, imputation_window_size: int = 3000000, k: int = 5,
                 threads: int = 4, verbose: bool = False):
        self.row_labels = row_labels
        self.sample_labels = sample_labels
        self.genomic_array = genomic_array
        assert isinstance(self.genomic_array, (np.ndarray, np.generic)), f'{self.genomic_array is None}'
        # check that matrix is properly oriented for imputation
        if genomic_array.shape[0] < genomic_array.shape[1]:
            raise RuntimeWarning(f'Number of Genomic Locations {genomic_array.shape[0]} < Number of Samples '
                                 f'{genomic_array.shape[1]}, Transpose Matrix if Rows Represent Samples')
        self.imputation_window_size = imputation_window_size
        self.k = k
        self.threads = threads
        self.verbose = verbose
        self.pool = None
        self.global_neighbors = None
        self.imputation_windows = None
        self.neighbor_dicts = None
        self.genomic_array_t = genomic_array.T

    def impute_windows(self):
        """Imputation run command"""
        self.get_global_neighbors()
        self.get_imputation_windows()
        self.get_window_values()
        self.watch_pool()
        self.set_missing_values()

    def get_global_neighbors(self):
        """Calculate global neighbors for sites without local information"""
        self.global_neighbors: List[List[float]] = get_euclidean(self.genomic_array_t)

    def get_imputation_windows(self):
        """Retrieve imputation windows based on input row labels"""
        self.imputation_windows = GenomeImputationWindows(site_labels=self.row_labels,
                                                          imputation_window_size=self.imputation_window_size)

    def get_window_values(self):
        """Impute pairwise distance for all imputation windows. Imputation windows are computed in parallel"""
        self.pool = mp.Pool(processes=self.threads)
        window_manager = mp.Manager()
        self.neighbor_dicts = window_manager.dict()
        for window_count, window in enumerate(self.imputation_windows.windows):
            if window:
                first_slice = self.row_labels.index(window[0])
                second_slice = self.row_labels.index(window[-1]) + 1
                window_values = np.copy(self.genomic_array[first_slice:second_slice, :])
                pairwise_distance_kwargs = dict(return_dict=self.neighbor_dicts,
                                                window_count=window_count,
                                                array=window_values,
                                                global_neighbors=self.global_neighbors,
                                                sample_labels=self.sample_labels,
                                                k=self.k)
                self.pool.apply_async(self.get_pairwise_distance, kwds=pairwise_distance_kwargs,
                                      error_callback=imputation_process_error)
            else:
                self.neighbor_dicts[window_count] = 'Empty Window'

    def watch_pool(self):
        """Watch pairwise distance operation and update status
        """
        completed_windows = None
        pbar = None
        if self.verbose:
            pbar = tqdm(total=len(self.imputation_windows.windows), desc='Calculating Window Pairwise Distance')
            completed_windows = len(self.neighbor_dicts)
            pbar.update(completed_windows)
        while len(self.neighbor_dicts) != len(self.imputation_windows.windows):
            if self.verbose:
                new_windows = len(self.neighbor_dicts) - completed_windows
                pbar.update(new_windows)
                completed_windows = len(self.neighbor_dicts)
        if self.verbose:
            new_windows = len(self.neighbor_dicts) - completed_windows
            pbar.update(new_windows)
            pbar.close()
        self.pool.close()

    @staticmethod
    def get_pairwise_distance(return_dict=None, window_count=None, global_neighbors=None,
                              array=None, sample_labels=None, k=None):
        """Method to calculate pairwise distance and return nearest neighbors. Distance matrix is discarded to
        reduce the memory overhead"""
        distance_matrix = get_euclidean(array.T, global_neighbors)
        local_neighbor_dict = {}
        for sample in zip(sample_labels, distance_matrix):
            sorted_sample_key = [a for b, a in sorted(zip(sample[1], sample_labels))]
            neighbors = [sample_labels.index(x) for x in sorted_sample_key[1:(k + 1)]]
            local_neighbor_dict[sample_labels.index(sample[0])] = neighbors
        return_dict[window_count] = local_neighbor_dict

    def set_missing_values(self):
        """Set imputed methylation values in genomic array and return results """
        # ensure all label are strings to set on label and not index
        imputed_values = []
        for label, row in tqdm(zip(self.row_labels, self.genomic_array), total=len(self.row_labels),
                               desc='Setting Imputed Values', disable=False if self.verbose else True):
            current_values = np.copy(row)
            null_methylation_values = np.argwhere(np.isnan(current_values))
            if len(null_methylation_values) > 0:
                null_site_neighbor_dict = self.neighbor_dicts[self.imputation_windows.site_window_dict[label]]
                for null_site in null_methylation_values:
                    neighbors = null_site_neighbor_dict[null_site[0]]
                    neighbor_average = round(np.nanmean(row[neighbors]), 5)
                    current_values[null_site[0]] = neighbor_average
            imputed_values.append(current_values)
        self.genomic_array = np.array(imputed_values)
