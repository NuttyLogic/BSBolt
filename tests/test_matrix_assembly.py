import random
import subprocess
import unittest

import numpy as np

from BSBolt.Matrix.MatrixAggregator import AggregateMatrix
from BSBolt.Utils.CGmapIterator import OpenCGmap
from BSBolt.Impute.Impute_Utils.ImputationFunctions import get_bsb_matrix
from tests.TestHelpers import test_directory


test_cgmap_file = [f'{test_directory}/TestData/test_cgmap_files/1.cgmap.gz',
                   f'{test_directory}/TestData/test_cgmap_files/2.cgmap',
                   f'{test_directory}/TestData/test_cgmap_files/3.cgmap',
                   f'{test_directory}/TestData/test_cgmap_files/4.cgmap']


def get_site_count(file_path):
    site_counts = {}
    for line in OpenCGmap(file_path):
        site_counts[f'{line[0]}:{line[2]}'] = int(line[7])
    return site_counts


all_counts = {}
for count, cgmap in enumerate(test_cgmap_file):
    all_counts[count] = get_site_count(cgmap)

low_count_sites = [key for key, count in all_counts[2].items() if count < 10]
missing_low_count = [site for site in low_count_sites if site not in all_counts[1]]

val_sites = set(random.sample(missing_low_count, 100))

# missing and low count sites should be excluded in matrix
missing_and_low_count = AggregateMatrix(file_list=test_cgmap_file[0:3], min_site_coverage=10,
                                        site_proportion_threshold=.9, verbose=False, threads=4)
missing_and_low_count.aggregate_matrix()

# assess count matrix
missing_and_low_count_count = AggregateMatrix(file_list=test_cgmap_file[0:3], min_site_coverage=10,
                                              site_proportion_threshold=.9, verbose=False, threads=4,
                                              count_matrix=True)
missing_and_low_count_count.aggregate_matrix()
# extract methylation values from count matrix
meth_count = missing_and_low_count_count.meth_matrix[:, 0::2]
total_counts = missing_and_low_count_count.meth_matrix[:, 1::2]
missing_and_low_count_meth_values = meth_count / total_counts


# missing sites with low count should be included in matrix
missing_and_five_count = AggregateMatrix(file_list=test_cgmap_file[0:3], min_site_coverage=5,
                                         site_proportion_threshold=.9, verbose=False, threads=4)
missing_and_five_count.aggregate_matrix()

# cg only sites
cg_only_test = AggregateMatrix(file_list=test_cgmap_file, min_site_coverage=5,
                               site_proportion_threshold=.9, verbose=False, cg_only=True, threads=4)
cg_only_test.aggregate_matrix()


# cg only high proportion
cg_only_test_high = AggregateMatrix(file_list=test_cgmap_file, min_site_coverage=10,
                                    site_proportion_threshold=1, verbose=False, cg_only=True, threads=4)
cg_only_test_high.aggregate_matrix()

# test_command line
test_matrix_output = f'{test_directory}/TestData/test_cgmap_files/matrix_test.txt'

bsb_matrix_commands = ['python3', '-m', 'BSBolt', 'AggregateMatrix',
                       '-F', f'{test_cgmap_file[0]},{test_cgmap_file[1]},{test_cgmap_file[2]}',
                       '-S', f'S1,S2,S3', '-O', test_matrix_output,
                       '-verbose', '-min-coverage', '10', '-min-sample', '0.9', '-t', '4']
subprocess.run(bsb_matrix_commands)

test_matrix, test_site_order, test_samples = get_bsb_matrix(test_matrix_output)

bsb_count_matrix_commands = ['python3', '-m', 'BSBolt', 'AggregateMatrix',
                             '-F', f'{test_cgmap_file[0]},{test_cgmap_file[1]},{test_cgmap_file[2]}',
                             '-S', f'S1,S2,S3', '-O', f'{test_matrix_output}_count', '-count', '-min-coverage',
                             '10', '-min-sample', '0.9', '-t', '4']

subprocess.run(bsb_count_matrix_commands)


class TestMatrixAggregation(unittest.TestCase):

    def setUp(self):
        pass

    def test_site_missing_low_count(self):
        # test site missing in 2.cgmap and low coverage in 3.cgmap are not included
        for site in val_sites:
            self.assertNotIn(site, missing_and_low_count.matrix_sites)

    def test_value_to_count_matrix(self):
        for value_values, count_values in zip(missing_and_low_count.meth_matrix, missing_and_low_count_meth_values):
            if np.isnan(sum(value_values)):
                # assert both sums are nan
                self.assertEqual(sum([np.isnan(sum(value_values)), np.isnan(sum(count_values))]), 2)
            else:
                self.assertAlmostEqual(sum(value_values), sum(count_values))

    def test_five_count(self):
        # test sites missing in 2.cgmap and low coverage in 3.cgmap are included with lower read count threshold
        for site in val_sites:
            self.assertIn(site, missing_and_five_count.matrix_sites)

    def test_cg_only(self):
        # test only cg sites are included
        self.assertEqual(len(cg_only_test.matrix_sites), 18842)

    def test_cg_only_high(self):
        # test missing and low count cg sites are not included in final matrix
        self.assertEqual(len(cg_only_test_high.matrix_sites), 16170)

    def test_command_line(self):
        for site in val_sites:
            self.assertNotIn(site, test_site_order)


if __name__ == '__main__':
    unittest.main()
