import gzip
import io
from typing import List, TextIO, Tuple, Union
import numpy as np
from bsbolt.GenotypeMatrix.GenotypeSiteCounter import GenotypeSiteCollector
from bsbolt.GenotypeMatrix.GenotypeSiteAggregator import GenotypeSiteAggregator


def propagate_error(error):
    raise error


class GenotypeAggregateMatrix:
    """
    Aggregate Variant BED files into combined variant matrix. Variant files are first iterated through to count sites and
    then iterated to through to retrieve values. While slower this prevents the construction of large sparse matrices.
    Assembly is multi-threaded to improve performance.

    Params:

   * *file_list (list)*: list of file Variant BED files
   * *sample_list (list)*: if passed, sample labels for Variant BED files,
                           else labels taken from sample names
   * *min_site_coverage (int)*: minimum read coverage for a CpG site to be considered for matrix, [10]
   * *site_proportion_threshold (float)*: proportion of samples that must have valid non-null
                                         variant calls for a site to be included in matrix, [0.9]
   * *output_path (str)*: path to output file
   * *verbose (bool)*: verbose matrix assembly, [False]
   * *threads (int)*: threads available for matrix aggregation, [1]
   * *encoding (bool)*: encoded matrix output, [False]

    Attributes:

    * *self.sample_list (list)*: list of samples as they are ordered in the variation matrix
    * *self.matrix_sites (tuple)*: tuple of ordered sites appearing in variant matrix, only set if no output path
                                   is provided
    * *self.gene_matrix (np.array)*: array of variant values (rows) by sample (columns), only set if not output path
                                    is provided

    Usage:

    ```python
    matrix = AggregateMatrix(**kwargs)
    matrix.aggregate_matrix()

    # access samples
    matrix.sample_list
    # access site list
    matrix.matrix_sites
    # access gene matrix
    matrix.gene_matrix
    """

    def __init__(self, file_list: List[str] = None, sample_list: List[str] = None, min_site_log_pvalue: int = 10,
                 site_proportion_threshold: float = 0.9, output_path: str = None,
                 verbose: bool = True, threads: int = 1, encoding: bool = False):
        self.file_list = file_list
        self.sample_list = sample_list
        if not self.sample_list:
            self.sample_list = []
            for file in self.file_list:
                self.sample_list.append(file.replace('\n', '').split('/')[-1])
        self.min_log_pvalue = min_site_log_pvalue
        self.verbose = verbose
        self.site_proportion_threshold = site_proportion_threshold
        self.output_path = output_path
        self.threads = threads
        self.gene_matrix = None
        self.matrix_sites = None
        self.encoding = encoding

    def aggregate_matrix(self):
        """Iterate through passed Variant files and aggregate matrix"""
        matrix_sites = self.collect_matrix_sites()
        gene_matrix = self.assemble_matrix(matrix_sites)
        if self.output_path:
            self.output_matrix(gene_matrix, matrix_sites)
        else:
            self.gene_matrix = gene_matrix
            self.matrix_sites = matrix_sites

    def collect_matrix_sites(self) -> Tuple[str]:
        """Iterate through individual files to get consensus site counts"""
        site_collector = GenotypeSiteCollector(variant_files=self.file_list,
                                            min_site_log_pvalue=self.min_log_pvalue,
                                            verbose=self.verbose,
                                            threads=self.threads,
                                            encoding=self.encoding,
                                            site_proportion_threshold=int(self.site_proportion_threshold
                                                                          * len(self.sample_list)))

        matrix_sites = site_collector.collect_consensus_sites()
        matrix_sites.sort(key=lambda x: x.split(':')[1])
        matrix_sites.sort(key=lambda x: x.split(':')[0])
        return tuple(matrix_sites)

    def assemble_matrix(self, matrix_sites: List[str]) -> np.ndarray:
        """Append sites to site list"""
        site_aggregator = GenotypeSiteAggregator(variant_files=self.file_list,
                                              min_site_log_pvalue=self.min_log_pvalue,
                                              verbose=self.verbose,
                                              threads=self.threads,
                                              encoding = self.encoding)
        matrix, matrix_order = site_aggregator.assemble_matrix(matrix_sites=matrix_sites)
        self.sample_list = [self.sample_list[pos] for pos in matrix_order]
        return matrix

    def get_output_object(self) -> Union[TextIO, io.BufferedWriter]:
        if self.output_path.endswith('gz'):
            return io.BufferedWriter(gzip.open(self.output_path, 'wb'))
        return open(self.output_path, 'w')

    def output_matrix(self, gene_matrix: np.ndarray, matrix_sites: List[str]):
        """Output sorted aggregated matrix"""
        out = self.get_output_object()
        with out as matrix:
            sample_labels = '\t'.join([str(sample) for sample in self.sample_list])
            matrix.write(f'Site\t{sample_labels}\n')
            for site, site_values in zip(matrix_sites, gene_matrix):
                gene_values = '\t'.join([f'{value}' for value in site_values])
                matrix.write(f'{site}\t{gene_values}\n')
