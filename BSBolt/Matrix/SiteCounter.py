import multiprocessing as mp
from typing import Deque, List
from BSBolt.Utils.CGmapIterator import OpenCGmap
from BSBolt.Utils.UtilityFunctions import propagate_error
from tqdm import tqdm


class CGmapSiteCollector:
    """Iterate through provided CGmap files counting the number of times a specific site is observed between all the
       cgmap files.
       Keyword Arguments
           cgmap_files (list): list of cgmap file paths
           min_site_coverage (int): min site coverage for a valid observation
           cg_only (bool): only consider cg context sites
           batch_size (int): chunk size for methylation observations to return
           verbose (bool): tqdm output
           threads (int): number of threads available
    """

    def __init__(self, cgmap_files: List[str] = None,
                 min_site_coverage: int = 1,
                 cg_only: bool = False,
                 batch_size: int = 40000,
                 verbose: bool = False,
                 threads: int = 8):
        self.cgmap_files = cgmap_files
        self.min_site_coverage = min_site_coverage
        self.cg_only = cg_only
        self.batch_size = batch_size
        self.verbose = verbose
        self.threads = threads
        self.observed_site_count = {}

    def collect_consensus_sites(self):
        # initialize collection_pool and manager queue
        collection_pool = mp.Pool(processes=self.threads)
        manager = mp.Manager()
        collection_queue = manager.Queue(maxsize=50)
        completed_samples = manager.list()
        for cgmap in self.cgmap_files:
            collection_pool.apply_async(self.collect_methylation_sites,
                                        kwds=dict(cgmap_path=cgmap, return_queue=collection_queue,
                                                  completed_samples=completed_samples),
                                        error_callback=propagate_error)
        collection_pool.close()
        # if verbose initialize tqdm progress bar
        update_count = 0
        pbar = tqdm(total=len(self.cgmap_files), desc='Collecting Methylation Sites') if self.verbose else None
        while True:
            if len(completed_samples) == len(self.cgmap_files) and collection_queue.empty():
                break
            sites = collection_queue.get(block=True)
            for site in sites:
                if site in self.observed_site_count:
                    self.observed_site_count[site] += 1
                else:
                    self.observed_site_count[site] = 1
            if pbar:
                update_number = len(completed_samples) - update_count
                pbar.update(update_number)
                update_count = len(completed_samples)
        if pbar:
            update_number = len(completed_samples) - update_count
            pbar.update(update_number)
            pbar.close()

    def process_sites(self, sites: list):
        for site in sites:
            if site in self.observed_site_count:
                self.observed_site_count[site] += 1
            else:
                self.observed_site_count[site] += 1

    def collect_methylation_sites(self, cgmap_path: str = None, return_queue: Deque = None,
                                  completed_samples: List = None):
        site_count, methylation_sites = 0, []
        for line_info in OpenCGmap(cgmap_path):
            try:
                line_info[7]
            except IndexError:
                break
            else:
                if self.check_cg(line_info[4]):
                    if int(line_info[7]) >= self.min_site_coverage:
                        methylation_sites.append(f'{line_info[0]}:{line_info[2]}')
                        site_count += 1
            if site_count == 10000:
                return_queue.put(tuple(methylation_sites), block=True)
                site_count, methylation_sites = 0, []
        if methylation_sites:
            return_queue.put(tuple(methylation_sites), block=True)
        completed_samples.append('done')

    def check_cg(self, nucleotide_context: str) -> bool:
        if self.cg_only:
            if nucleotide_context == 'CG':
                return True
            return False
        return True
