import os
import pickle
import queue
import time
import unittest

from BSBolt.CallMethylation.CallMethylationValues import CallMethylationValues
from BSBolt.CallMethylation.CallMethylationVector import CallMethylationVector
from tests.TestHelpers import bsb_directory, z_test_of_proportion


bsb_index = os.path.exists(f'{bsb_directory}tests/TestSimulations/BSB_pe.chr11.pkl')
bsb_alignment_file = os.path.exists(f'{bsb_directory}tests/BSB_pe_test.sorted.bam')

if not all([bsb_index, bsb_alignment_file]):
    print('Simulate date before running tests')
    raise ValueError

chr11_methylation_sites = {}

with open(f'{bsb_directory}tests/TestSimulations/BSB_pe.chr11.pkl', 'rb') as sim_sites:
    sim_values = pickle.load(sim_sites)
    chr11_methylation_sites.update(sim_values['Watson'])
    chr11_methylation_sites.update(sim_values['Crick'])


# call methylation values

value_queue = queue.Queue()

value_start = time.time()

value_caller = CallMethylationValues(input_file=f'{bsb_directory}tests/BSB_pe_test.sorted.bam',
                                     genome_database=f'{bsb_directory}tests/TestData/BSB_Test_DB/',
                                     contig='chr11', return_queue=value_queue, ignore_orphans=False,
                                     ignore_overlap=True, min_base_quality=1)
value_caller.call_methylation()

print(f'Values Called: {time.time() - value_start: .3f} sec')

all_values = []

while True:
    values = value_queue.get()
    all_values.extend(values)
    if value_queue.empty():
        break

# process methylation values

cleaned_values = {}

for value in all_values:
    site_label = f'{value[8]}:{value[7] - 1}'
    cleaned_values[site_label] = [value[1], value[2]]


vector_queue = queue.Queue()

vector_start = time.time()

vector_caller = CallMethylationVector(input_file=f'{bsb_directory}tests/BSB_pe_test.sorted.bam',
                                      genome_database=f'{bsb_directory}tests/TestData/BSB_Test_DB/',
                                      contig='chr11', return_queue=vector_queue, min_base_quality=1)
vector_caller.call_methylation()

print(f'Vectors Called: {time.time() - vector_start: .3f} sec')


all_vectors = []


while True:
    values = vector_queue.get()
    all_vectors.extend(values[1])
    if vector_queue.empty():
        break

vector_point_values = {}

vector_point_start = time.time()


def get_vector_point_values(vector, chromosome, value_dict):
    for meth_call, pos in zip(vector[3], vector[4]):
        if f'{chromosome}:{pos}' in value_dict:
            if meth_call:
                value_dict[f'{chromosome}:{pos}'][0] += 1
            else:
                value_dict[f'{chromosome}:{pos}'][1] += 1
        else:
            if meth_call:
                value_dict[f'{chromosome}:{pos}'] = [1, 0]
            else:
                value_dict[f'{chromosome}:{pos}'] = [0, 1]


# get point values
for vector in all_vectors:
    get_vector_point_values(vector, 'chr11', vector_point_values)

print(f'Vector Point Values Called: {time.time() - vector_point_start: .3f} sec')

vector_point_value_comparison = {}

mismatch_count = 0
differences = []
for site, vector_values in vector_point_values.items():
    pileup_values = cleaned_values[site]
    if vector_values != pileup_values:
        mismatch_count += 1
        difference = abs(vector_values[0] - pileup_values[0]) + abs(vector_values[1] - pileup_values[1])
        differences.append(difference)


site_comparisons = {}

for site, values in vector_point_values.items():
    site_comparison = dict(coverage_difference=0, simulation_beta=0, mapped_beta=0, beta_z_value=0)
    site_coverage = sum(values)
    reference_values = chr11_methylation_sites[site]
    ref_meth = int(reference_values[3])
    ref_unmeth = int(reference_values[4])
    reference_coverage = ref_meth + ref_unmeth
    try:
        site_comparison['coverage_difference'] = abs(site_coverage - reference_coverage)
        site_comparison['mapped_beta'] = values[0] / site_coverage
        site_comparison['simulation_beta'] = ref_meth / reference_coverage
        z = abs(z_test_of_proportion(a_yes=values[0], a_no=values[1], b_yes=ref_meth, b_no=ref_unmeth))
        site_comparison['beta_z_value'] = z
    except ZeroDivisionError:
        site_comparison[site] = dict(coverage_difference=0, simulation_beta=0, mapped_beta=0, beta_z_value=0)
    site_comparisons[site] = site_comparison


class TestMethylationCalling(unittest.TestCase):
    """
    Simulated sites will not always map uniquely to reference genome, set tight difference limits to test correct
    methylation calling, but not absolute differences
    """

    def setUp(self):
        pass

    def test_meth_point_vector(self):
        self.assertLess(len(differences), 5)

    def test_read_coverage(self):
        # set coverage difference between simulated and mapped reads to consider site out of tolerance
        coverage_difference_tolerance = 5
        # count number of sites out of tolerance
        out_of_tolerance_sites = 0
        for test_site in site_comparisons.values():
            if test_site['coverage_difference'] > coverage_difference_tolerance:
                out_of_tolerance_sites += 1
        self.assertLessEqual(out_of_tolerance_sites, 20)

    def test_beta_proportion(self):
        # set z threshold
        z_threshold = 3
        # count site with z score above threshold
        z_site_count = 0
        for test_site in site_comparisons.values():
            if test_site['beta_z_value'] >= z_threshold:
                z_site_count += 1
        self.assertLessEqual(z_site_count, 20)


if __name__ == '__main__':
    unittest.main()
