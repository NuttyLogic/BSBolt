import multiprocessing as mp
from random import sample
from typing import Dict, List
from bsbolt.Utils.CGmapIterator import OpenCGmap
from bsbolt.Utils.UtilityFunctions import propagate_error
import numpy as np
from tqdm import tqdm


class GenotypeSiteAggregator:
    """Pull consensus variant sites into a variant matrix
       Keyword Arguments:
              variant_files (list): file paths
              min_site_coverage (int): minimum site coverage for a valid variant call
              verbose (bool): tqdm output
              threads (int): threads available for aggregation
              encoding (bool): output encoded matrix
        """

    def __init__(self, variant_files: List[str] = None,
                 min_site_log_pvalue: int = 1,
                 verbose: bool = False,
                 encoding: bool = False,
                 threads: int = 8):
        self.variant_files = variant_files
        self.min_site_log_pvalue = min_site_log_pvalue
        self.verbose = verbose
        self.threads = threads
        self.encoding = encoding

    def assemble_matrix(self, matrix_sites: List[str] = None):
        # initialize collection_pool and manager queue
        matrix_pool = mp.Pool(processes=self.threads)
        manager = mp.Manager()
        matrix_queue = manager.Queue(maxsize=10)
        completed_samples = manager.list()
        # initialize empty array of zeros, fill with null values for safety
        for count, variant in enumerate(self.variant_files):
            matrix_pool.apply_async(self.collect_variant_sites,
                                    kwds=dict(sample_pos=count, variant_path=variant, return_queue=matrix_queue,
                                              completed_samples=completed_samples, matrix_sites=matrix_sites),
                                    error_callback=propagate_error)
        matrix_pool.close()
        # if verbose initialize tqdm progress bar
        update_count = 0
        pbar = tqdm(total=len(self.variant_files), desc='Setting Variant Sites') if self.verbose else None
        variant_matrix, matrix_order = [], []
        while True:
            if len(completed_samples) == len(self.variant_files) and matrix_queue.empty():
                break
            sample_pos, sample_values = matrix_queue.get(block=True)
            matrix_order.append(sample_pos)
            if len(sample_values.shape) > 1:
                variant_matrix.append(sample_values[:, 0])
                variant_matrix.append(sample_values[:, 1])
            else:
                variant_matrix.append(sample_values)
            if pbar:
                update_number = len(completed_samples) - update_count
                pbar.update(update_number)
                update_count = len(completed_samples)
        if pbar:
            update_number = len(completed_samples) - update_count
            pbar.update(update_number)
            pbar.close()
        return np.array(variant_matrix).T, matrix_order

    def collect_variant_sites(self, sample_pos: int = None, variant_path: str = None, return_queue: mp.Queue = None,
                                  completed_samples: List = None, matrix_sites: Dict = None):
        # initialize an array for all sites and fill with null values
        matrix_key = {site: count for count, site in enumerate(matrix_sites)}
        sample_array = np.zeros((len(matrix_sites)),'U3')
        sample_array.fill(np.nan)
        for line_info in OpenCGmap(variant_path):
            try:
                line_info[4]
            except IndexError:
                break
            else:
                if float(line_info[4]) >= self.min_site_log_pvalue:
                    site_label = f'{line_info[0]}:{line_info[1]}:{line_info[5]}'
                    if site_label in matrix_key:
                        if self.encoding:
                            ref = line_info[5]
                            if line_info[6] == ref:
                                line_info[6] = 0
                            elif "," not in line_info[6]:
                                line_info[6] = 2
                            else:
                                line_info[6] = 1
                        sample_array[matrix_key[site_label]] = line_info[6]
        return_queue.put((sample_pos, sample_array), block=True)
        completed_samples.append('done')
