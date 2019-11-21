#! /usr/env python3

import gzip
import io
from BSBolt.Matrix.SiteCounter import CGmapSiteCollector
from BSBolt.Matrix.SiteAggregator import CGmapSiteAggregator


def propagate_error(error):
    raise error


class AggregateMatrix:
    """
    Class to aggregate CGMap files into combined methylation matrix
    Keyword Arguments:
        file_list (list): list of file CGmap files
        sample_list (list): if passed, sample labels for CGmaps files, 
                            else labels taken from sample names
        min_site_coverage (int): minimum read coverage for a CpG site to be considered for matrix
        site_proportion_threshold (float): proportion of samples that must have valid non-null
                                           methylation calls for a site to be included in matrix
        output_path (str): path to output file
        cg_only (bool): consider all cystosines or only CpG sites
        verbose (bool): verbose matrix assembly, tqdm output
        threads (int): number of threads available for matrix aggregation
    Attributes:
       self.file_list (list): list of file CGmap files
       self.sample_list (list): if passed, sample labels for CGmaps files, 
                           else labels taken from sample names
       self.min_site_coverage (int): minimum read coverage for a CpG site to be considered for matrix
       self.site_proportion_threshold (float): proportion of samples that must have valid non-null
                                          methylation calls for a site to be included in matrix
       self.output_path (str): path to output file
       self.cg_only (bool): consider all cystosines or only CpG sites
       self.disable_tqdm (bool): disable tqdm
       self.threads (int): number of threads available for matrix aggregation
       self.matrix_sites (tuple): tuple of ordered sites appearing in methylation matrix, only set if no output path
                                  is provided
       self.meth_matrix (np.array): array of methylation values (rows) by sample (columns), only set if not output path
                                    is provided
    """

    def __init__(self, file_list=None, sample_list=None, min_site_coverage=10,
                 site_proportion_threshold=0.9, output_path=None, cg_only=False, verbose=True, threads=1):
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
        self.meth_matrix = None
        self.matrix_sites = None

    def aggregate_matrix(self):
        """Collapse matrix"""
        matrix_sites = self.collect_matrix_sites()
        meth_matrix = self.assemble_matrix(matrix_sites)
        if self.output_path:
            self.output_matrix(meth_matrix, matrix_sites)
        else:
            self.meth_matrix = meth_matrix
            self.matrix_sites = matrix_sites

    @staticmethod
    def sort_sites(key_list):
        """Takes list of chrx:XXXX sites and sorts them based on chr then genomic location"""
        sorting_list = [[], []]
        for site in key_list:
            site_split = site.split(':')
            sorting_list[0].append(site_split[0])
            sorting_list[1].append(int(site_split[1]))
        key_list = [x for (z, y, x) in sorted(zip(sorting_list[0], sorting_list[1], key_list),
                                              key=lambda sorting_values: (sorting_values[0], sorting_values[1]))]
        return key_list

    def collect_matrix_sites(self):
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

    def assemble_matrix(self, matrix_sites):
        """Append sites to site list"""
        site_aggregator = CGmapSiteAggregator(cgmap_files=self.file_list,
                                              min_site_coverage=self.min_coverage,
                                              verbose=self.verbose,
                                              threads=self.threads)
        site_aggregator.assemble_matrix(matrix_sites=matrix_sites)
        return site_aggregator.meth_matrix

    def get_output_object(self):
        if self.output_path.endswith('gz'):
            return io.BufferedWriter(gzip.open(self.output_path, 'wb'))
        return open(self.output_path, 'w')

    def output_matrix(self, meth_matrix, matrix_sites):
        """Output sorted aggregated matrix"""
        out = self.get_output_object()
        with out as matrix:
            sample_labels = '\t'.join([str(sample) for sample in self.sample_list])
            matrix.write(f'Site\t{sample_labels}\n')
            for site, site_values in zip(matrix_sites, meth_matrix):
                meth_values = '\t'.join([f'{value:.4f}' for value in site_values])
                matrix.write(f'{site}\t{meth_values}\n')
