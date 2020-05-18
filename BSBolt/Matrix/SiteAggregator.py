import multiprocessing as mp
from typing import Dict, Deque, List
from BSBolt.Utils.CGmapIterator import OpenCGmap
from BSBolt.Utils.UtilityFunctions import propagate_error
import numpy as np
from tqdm import tqdm


class CGmapSiteAggregator:
    """Pull consensus methylation sites into a methylation matrix
       Keyword Arguments:
              cgmap_files (list): file paths
              min_site_coverage (int): minimum site coverage for a valid methylation call
              verbose (bool): tqdm output
              threads (int): threads available for aggregation
        """

    def __init__(self, cgmap_files: List[str] = None,
                 min_site_coverage: int = 1,
                 verbose: bool = False,
                 threads: int = 8,
                 count_matrix: bool = False):
        self.cgmap_files = cgmap_files
        self.min_site_coverage = min_site_coverage
        self.verbose = verbose
        self.threads = threads
        self.count_matrix = count_matrix
        self.meth_matrix = None

    def assemble_matrix(self, matrix_sites: Dict = None):
        # initialize collection_pool and manager queue
        matrix_pool = mp.Pool(processes=self.threads)
        manager = mp.Manager()
        matrix_queue = manager.Queue(maxsize=10)
        completed_samples = manager.list()
        # initialize empty array of zeros, fill with null values for safety
        if self.count_matrix:
            self.meth_matrix = np.zeros((len(matrix_sites), len(self.cgmap_files) * 2))
        else:
            self.meth_matrix = np.zeros((len(matrix_sites), len(self.cgmap_files)))
        for count, cgmap in enumerate(self.cgmap_files):
            matrix_pool.apply_async(self.collect_methylation_sites,
                                    kwds=dict(sample_pos=count, cgmap_path=cgmap, return_queue=matrix_queue,
                                              completed_samples=completed_samples, matrix_sites=matrix_sites),
                                    error_callback=propagate_error)
        matrix_pool.close()
        # if verbose initialize tqdm progress bar
        update_count = 0
        pbar = tqdm(total=len(self.cgmap_files), desc='Setting Methylation Sites') if self.verbose else None
        while True:
            if len(completed_samples) == len(self.cgmap_files) and matrix_queue.empty():
                break
            sample_pos, sample_values = matrix_queue.get(block=True)
            sample_span = sample_pos
            if self.count_matrix:
                sample_span = (sample_pos * 2, sample_pos * 2 + 1)
            self.meth_matrix[:, sample_span] = sample_values
            if pbar:
                update_number = len(completed_samples) - update_count
                pbar.update(update_number)
                update_count = len(completed_samples)
        if pbar:
            update_number = len(completed_samples) - update_count
            pbar.update(update_number)
            pbar.close()

    def collect_methylation_sites(self, sample_pos: int = None, cgmap_path: str = None, return_queue: Deque = None,
                                  completed_samples: List = None, matrix_sites: Dict = None):
        # initialize an array for all sites and fill with null values
        if self.count_matrix:
            sample_array = np.zeros((len(matrix_sites), 2))
        else:
            sample_array = np.zeros((len(matrix_sites)))
        sample_array.fill(np.nan)
        for line_info in OpenCGmap(cgmap_path):
            try:
                line_info[7]
            except IndexError:
                break
            else:
                if int(line_info[7]) >= self.min_site_coverage:
                    site_label = f'{line_info[0]}:{line_info[2]}'
                    if site_label in matrix_sites:
                        if self.count_matrix:
                            sample_array[matrix_sites[site_label], :] = np.array((int(line_info[6]), int(line_info[7])))
                        else:
                            sample_array[matrix_sites[site_label]] = float(line_info[5])
        return_queue.put((sample_pos, sample_array), block=True)
        completed_samples.append('done')
