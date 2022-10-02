from collections import Counter
import multiprocessing as mp
from typing import List
from bsbolt.Utils.CGmapIterator import OpenCGmap
from bsbolt.Utils.UtilityFunctions import propagate_error
from tqdm import tqdm


class GenotypeSiteCollector:
    """Iterate through provided variant files counting the number of times a specific site is observed between all the
       variant files.
       Keyword Arguments
           variant_files (list): list of variant BED file paths
           min_site_coverage (int): min site coverage for a valid observation
           batch_size (int): chunk size for variant observations to return
           verbose (bool): tqdm output
           threads (int): number of threads available
           site_proportion_threshold (int): number of times a site must be observed to be included in matrix
    """

    def __init__(self, variant_files: List[str] = None,
                 min_site_log_pvalue: int = 1,
                 batch_size: int = 10000,
                 verbose: bool = False,
                 threads: int = 8,
                 encoding: bool = False,
                 site_proportion_threshold: int = 1):
        self.variant_files = variant_files
        self.min_site_log_pvalue = min_site_log_pvalue
        self.batch_size = batch_size
        self.verbose = verbose
        self.threads = threads
        self.encoding = encoding
        self.site_proportion_threshold = site_proportion_threshold
        self.observed_site_count = Counter()

    def collect_consensus_sites(self) -> List[str]:
        # initialize collection_pool and manager queue
        collection_pool = mp.Pool(processes=self.threads)
        manager = mp.Manager()
        collection_queue = manager.Queue(maxsize=10)
        completed_samples = manager.list()
        for variant in self.variant_files:
            collection_pool.apply_async(self.collect_variant_sites,
                                        kwds=dict(variant_path=variant, return_queue=collection_queue,
                                                  completed_samples=completed_samples),
                                        error_callback=propagate_error)
        collection_pool.close()
        # if verbose initialize tqdm progress bar
        update_count = 0
        pbar = tqdm(total=len(self.variant_files), desc='Collecting Variant Sites') if self.verbose else None
        while True:
            if len(completed_samples) == len(self.variant_files) and collection_queue.empty():
                break
            sites = collection_queue.get(block=True)
            self.observed_site_count.update(sites)
            if pbar:
                update_number = len(completed_samples) - update_count
                pbar.update(update_number)
                update_count = len(completed_samples)
        if pbar:
            update_number = len(completed_samples) - update_count
            pbar.update(update_number)
            pbar.close()
        print(self.site_proportion_threshold)
        return [site for site, count in self.observed_site_count.items() if count >= self.site_proportion_threshold]

    def collect_variant_sites(self, variant_path: str = None, return_queue: mp.Queue = None,
                                  completed_samples: List = None):
        site_count, variant_sites = 0, []
        for line_info in OpenCGmap(variant_path):
            try:
                line_info[4]
            except IndexError:
                break
            else:
                if float(line_info[4]) >= self.min_site_log_pvalue:
                    variant_sites.append(f'{line_info[0]}:{line_info[1]}:{line_info[5]}')
                    site_count += 1
            if site_count == self.batch_size:
                return_queue.put(tuple(variant_sites), block=True)
                site_count, variant_sites = 0, []
        if variant_sites:
            return_queue.put(tuple(variant_sites), block=True)
        completed_samples.append('done')

