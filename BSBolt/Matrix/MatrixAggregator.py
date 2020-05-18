import gzip
import io
from typing import Dict, List, TextIO, Union
import numpy as np
from BSBolt.Matrix.SiteCounter import CGmapSiteCollector
from BSBolt.Matrix.SiteAggregator import CGmapSiteAggregator


def propagate_error(error):
    raise error


class AggregateMatrix:
    """
    Aggregate CGMap files into combined methylation matrix. CGmap files are first iterated through to count sites and
    then iterated to through to retrieve values. While slower this prevents the construction of large sparse matrices.
    Assembly is multi-threaded to improve performance.

    Params:

   * *file_list (list)*: list of file CGmap files
   * *sample_list (list)*: if passed, sample labels for CGmaps files,
                           else labels taken from sample names
   * *min_site_coverage (int)*: minimum read coverage for a CpG site to be considered for matrix, [10]
   * *site_proportion_threshold (float)*: proportion of samples that must have valid non-null
                                         methylation calls for a site to be included in matrix, [0.9]
   * *output_path (str)*: path to output file
   * *cg_only (bool)*: consider all cytosines or only CpG sites, [False]
   * *verbose (bool)*: verbose matrix assembly, [False]
   * *threads (int)*: threads available for matrix aggregation, [1]

    Attributes:

    * *self.sample_list (list)*: list of samples as they are ordered in the methylation matrix
    * *self.matrix_sites (tuple)*: tuple of ordered sites appearing in methylation matrix, only set if no output path
                                   is provided
    * *self.meth_matrix (np.array)*: array of methylation values (rows) by sample (columns), only set if not output path
                                    is provided

    Usage:

    ```python
    matrix = AggregateMatrix(**kwargs)
    matrix.aggregate_matrix()

    # access samples
    matrix.sample_list
    # access site list
    matrix.matrix_sites
    # access methylation matrix
    matrix.meth_matrix
    """

    def __init__(self, file_list: List[str] = None, sample_list: List[str] = None, min_site_coverage: int = 10,
                 site_proportion_threshold: float = 0.9, output_path: str = None, cg_only: bool = False,
                 verbose: bool = True, threads: int = 1, count_matrix: bool = False):
        self.file_list = file_list
        self.sample_list = sample_list
        if not self.sample_list:
            self.sample_list = []
            for file in self.file_list:
                self.sample_list.append(file.replace('\n', '').split('/')[-1])
        self.min_coverage = min_site_coverage
        self.verbose = verbose
        self.site_proportion_threshold = site_proportion_threshold
        self.output_path = output_path
        self.cg_only = cg_only
        self.threads = threads
        self.count_matrix = count_matrix
        self.meth_matrix = None
        self.matrix_sites = None

    def aggregate_matrix(self):
        """Iterate through passed CGmap files and aggregate matrix"""
        matrix_sites = self.collect_matrix_sites()
        meth_matrix = self.assemble_matrix(matrix_sites)
        if self.output_path:
            self.output_matrix(meth_matrix, matrix_sites)
        else:
            self.meth_matrix = meth_matrix
            self.matrix_sites = matrix_sites

    def collect_matrix_sites(self) -> Dict[str, int]:
        """Iterate through individual files to get consensus site counts"""
        site_collector = CGmapSiteCollector(cgmap_files=self.file_list,
                                            min_site_coverage=self.min_coverage,
                                            cg_only=self.cg_only,
                                            verbose=self.verbose,
                                            threads=self.threads)
        site_collector.collect_consensus_sites()
        matrix_sites = []
        site_count_threshold = int(self.site_proportion_threshold * len(self.sample_list))
        for site, count in site_collector.observed_site_count.items():
            if count >= site_count_threshold:
                matrix_sites.append(site)
        matrix_sites.sort(key=lambda x: int(x.split(':')[1]))
        matrix_sites.sort(key=lambda x: x.split(':')[0])
        return {site: count for count, site in enumerate(matrix_sites)}

    def assemble_matrix(self, matrix_sites: Dict[str, int]) -> np.ndarray:
        """Append sites to site list"""
        site_aggregator = CGmapSiteAggregator(cgmap_files=self.file_list,
                                              min_site_coverage=self.min_coverage,
                                              verbose=self.verbose,
                                              threads=self.threads,
                                              count_matrix=self.count_matrix)
        site_aggregator.assemble_matrix(matrix_sites=matrix_sites)
        return site_aggregator.meth_matrix

    def get_output_object(self) -> Union[TextIO, io.BufferedWriter]:
        if self.output_path.endswith('gz'):
            return io.BufferedWriter(gzip.open(self.output_path, 'wb'))
        return open(self.output_path, 'w')

    def output_matrix(self, meth_matrix: np.ndarray, matrix_sites: Dict[str, int]):
        """Output sorted aggregated matrix"""
        out = self.get_output_object()
        with out as matrix:
            sample_labels = '\t'.join([str(sample) for sample in self.sample_list])
            if self.count_matrix:
                sample_labels = '\t'.join([f'{sample}_meth_cytosine\t{sample}_total_cytosine'
                                          for sample in self.sample_list])
            matrix.write(f'Site\t{sample_labels}\n')
            for site, site_values in zip(matrix_sites, meth_matrix):
                if self.count_matrix:
                    meth_values = '\t'.join([f'{value}' for value in site_values])
                else:
                    meth_values = '\t'.join([f'{value:.6f}' for value in site_values])
                matrix.write(f'{site}\t{meth_values}\n')
