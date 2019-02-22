#! /usr/env python3

import numpy as np
from tqdm import tqdm
from BSB.BSB_Utils.CGmapIterator import OpenCGmap


class AggregateMatrix:
    """
    Class to aggregate CGMap files into combined methylation matrix
    Keyword Arguments:
        file_list (list): list of file CGmap files
        sample_list (list): if passed, sample labels for CGmaps files, 
                            else labels taken from sample names
        min_site_coverage (int): minimum read coverage for a CpG site to be considered for matrix
        site_proportion_threshold (float): proprotion of samples that must have valid non-null 
                                           methylation calls for a site to be included in matrix
        output_path (str): path to output file
        cg_only (bool): consider all cystosines or only CpG sites
        verbose (bool): verbose matrix assembly, tqdm output
    Attributes:
       self.file_list (list): list of file CGmap files
       self.sample_list (list): if passed, sample labels for CGmaps files, 
                           else labels taken from sample names
       self.min_site_coverage (int): minimum read coverage for a CpG site to be considered for matrix
       self.site_proportion_threshold (float): proprotion of samples that must have valid non-null 
                                          methylation calls for a site to be included in matrix
       self.output_path (str): path to output file
       self.cg_only (bool): consider all cystosines or only CpG sites
       self.disable_tqdm (bool): disable tqdm
       self.site_dict (dict): dict to store methylation values 
    """

    def __init__(self, file_list=None, sample_list=None, min_site_coverage=10,
                 site_proportion_threshold=0.9, output_path=None, cg_only=True, verbose=True):
        self.file_list = file_list
        self.sample_list = sample_list
        if not self.sample_list:
            self.sample_list = []
            for file in self.file_list:
                self.sample_list.append(file.split('/')[-1])
        self.min_coverage = min_site_coverage
        self.site_dict = None
        self.collapsed_matrix = None
        self.disable_tqdm = False if verbose else True
        self.site_proportion_threshold = site_proportion_threshold
        self.output_path = output_path
        self.cg_only = cg_only
        self.site_dict = {}

    def aggregate_matrix(self):
        """Collapse matrix"""
        self.get_site_counts()
        self.get_matrix_sites()
        self.set_matrix_sites()
        if self.output_path:
            self.output_matrix()

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

    def get_site_counts(self):
        """Iterate through individual files to get consensus site counts"""
        for sample_number, file in tqdm(enumerate(self.file_list), desc='Counting CpG Sites',
                                        total=len(self.file_list)):
            for line_info in OpenCGmap(cgmap=file):
                try:
                    line_info[7]
                except IndexError:
                    break
                if self.check_cg(line_info[4]):
                    if int(line_info[7]) >= self.min_coverage:
                        site_label = f'{line_info[0]}:{line_info[2]}'
                        try:
                            self.site_dict[site_label] += 1
                        except KeyError:
                            self.site_dict[site_label] = 1

    def get_matrix_sites(self):
        """Set sites seen across samples sbove proportion """
        site_count_threshold = int(self.site_proportion_threshold * len(self.sample_list))
        self.collapsed_matrix = {}
        for key, value in tqdm(self.site_dict.items(), desc='Getting Methylation Sites', total=len(self.site_dict)):
            if value >= site_count_threshold:
                site_values: np.array = np.zeros(len(self.sample_list))
                site_values.fill(np.nan)
                self.collapsed_matrix[key] = site_values
        del self.site_dict

    def set_matrix_sites(self):
        """Append sites to site list"""
        for sample_index, file in tqdm(enumerate(self.file_list), desc='Setting Methylation Values',
                                       total=len(self.file_list)):
            for line_info in OpenCGmap(cgmap=file):
                site_label = f'{line_info[0]}:{line_info[2]}'
                if int(line_info[7]) >= self.min_coverage and site_label in self.collapsed_matrix:
                    self.collapsed_matrix[site_label][sample_index] = float(line_info[5])

    def check_cg(self, nucleotide_context):
        if self.check_cg:
            if nucleotide_context == 'CG':
                return True
            return False
        return True

    def output_matrix(self):
        """Output sorted aggregated matrix"""
        print('Writing Matrix')
        key_list = self.sort_sites(list(self.collapsed_matrix.keys()))
        out = open(self.output_path, 'w')
        out.write(' \t%s\n' % '\t'.join(self.sample_list))
        for site in key_list:
            meth_values = self.collapsed_matrix[site]
            out.write('%s\t%s\n' % (site, '\t'.join(meth_values)))
        out.close()
