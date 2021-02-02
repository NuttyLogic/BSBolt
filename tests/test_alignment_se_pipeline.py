import gzip
import os
import pickle
import subprocess
import unittest

from bsbolt.Utils.AlignmentEvaluation import AlignmentEvaluator
from tests.TestHelpers import bsb_directory, z_test_of_proportion

# generate simulated reads
bsb_simulate_commands = ['python3', '-m', 'bsbolt', 'Simulate',
                         '-G', f'{bsb_directory}tests/TestData/BSB_test.fa', '-U',
                         '-O', f'{bsb_directory}tests/TestSimulations/BSB_se',
                         '-MR', '0.01', '-verbose', '-overwrite', '-RD', '20']
subprocess.run(bsb_simulate_commands)
print('Reads Simulated')
# map simulated reads

print('Building Methylation Index')
bsb_index_commands = ['python3', '-m', 'bsbolt', 'Index', '-G', f'{bsb_directory}tests/TestData/BSB_test.fa',
                      '-DB', f'{bsb_directory}tests/TestData/BSB_Test_DB']
subprocess.run(bsb_index_commands)
print('bsbolt Index Built')


bsb_align_commands = ['python3', '-m', 'bsbolt', 'Align',
                      '-DB', f'{bsb_directory}tests/TestData/BSB_Test_DB', '-F1',
                      f'{bsb_directory}tests/TestSimulations/BSB_se_1.fq', '-O',
                      f'{bsb_directory}tests/BSB_se_test', '-t', '10', '-UN', '-SP', '0.1']

subprocess.run(bsb_align_commands)

sorted_output = f'{bsb_directory}tests/BSB_se_test.sorted.bam'
if os.path.exists(sorted_output):
    subprocess.run(['rm', sorted_output])
    subprocess.run(['rm', f'{sorted_output}.bai'])

subprocess.run(['python3', '-m', 'bsbolt', 'Sort', '-I', f'{bsb_directory}tests/BSB_se_test.bam',
                '-O', f'{bsb_directory}tests/BSB_se_test.sorted.bam'])

print('Calling Methylation')

bs_call_methylation_args = ['python3', '-m', 'bsbolt', 'CallMethylation', '-I',
                            f'{bsb_directory}tests/BSB_se_test.sorted.bam',
                            '-O', f'{bsb_directory}tests/BSB_se_test',
                            '-DB', f'{bsb_directory}tests/TestData/BSB_Test_DB',
                            '-t', '6', '-verbose', '-BQ', '10', '-MQ', '20']
subprocess.run(bs_call_methylation_args)
print('Methylation Values Called')

# retrieve reference and test alignments
test_alignments = f'{bsb_directory}tests/BSB_se_test.sorted.bam'
test_fastqs = [f'{bsb_directory}tests/TestSimulations/BSB_se_1.fq']

evaluator = AlignmentEvaluator(duplicated_regions={'chr10': (0, 5000), 'chr15': (0, 5000)},
                               matching_target_prop=.95)

print('Evaluating Alignment')
read_stats = evaluator.evaluate_alignment(test_alignments, test_fastqs)


# import methylation calling dict
print('Evaluating Methylation Calls')
all_methylation_sites = {}

output_chromosome = ['chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']
for chrom in output_chromosome:
    with open(f'{bsb_directory}tests/TestSimulations/BSB_se.{chrom}_values.pkl', 'rb') as sim_sites:
        chrom_sites = pickle.load(sim_sites)
        all_methylation_sites.update(chrom_sites)

# import CGmap Calls
cgmap_sites = {}

for line in gzip.open(f'{bsb_directory}tests/BSB_se_test.CGmap.gz', 'rb'):
    processed_line = line.decode('utf-8').replace('\n', '').split('\t')
    cgmap_sites[f'{processed_line[0]}:{int(processed_line[2]) - 1}'] = dict(nucleotide=processed_line[1],
                                                                            methylation_level=processed_line[5],
                                                                            context=processed_line[4],
                                                                            methylated_reads=processed_line[6],
                                                                            total_reads=processed_line[7])

# get line comparisons

site_comparisons = {}
missing_sites = []
for site, cgmap_values in cgmap_sites.items():
    site_comparison = dict(coverage_difference=0, simulation_beta=0, mapped_beta=0, beta_z_value=0)
    cgmap_site_coverage = int(cgmap_values['total_reads'])
    try:
        reference_values = all_methylation_sites[site]
    except KeyError:
        # don't reference methylation at variant sites
        missing_sites.append(site)
        continue
    reference_coverage = int(reference_values[0]) + int(reference_values[1])
    site_comparison['coverage_difference'] = abs(cgmap_site_coverage - reference_coverage)
    site_comparison['simulation_beta'] = reference_values[0] / reference_coverage
    site_comparison['mapped_beta'] = cgmap_values['methylation_level']
    cgmap_meth = int(cgmap_sites[site]['methylated_reads'])
    cgmap_unmeth = int(cgmap_sites[site]['total_reads']) - int(cgmap_sites[site]['methylated_reads'])
    ref_meth = int(reference_values[0])
    ref_unmeth = int(reference_values[1])
    try:
        z = abs(z_test_of_proportion(a_yes=cgmap_meth, a_no=cgmap_unmeth, b_yes=ref_meth, b_no=ref_unmeth))
    except ZeroDivisionError:
        z = 0
    site_comparison['beta_z_value'] = z
    site_comparisons[site] = site_comparison


class TestBSBPipeline(unittest.TestCase):
    """ The first 5000bp for chr10 are duplicated as chr15 in the simulation reference. These regions will have
    mixed methylation values and coverage values"""

    def setUp(self):
        pass

    def test_read_coverage(self):
        # set coverage difference between simulated and mapped reads to consider site out of tolerance
        coverage_difference_tolerance = 5
        # count number of sites out of tolerance
        out_of_tolerance_sites = 0
        for label, test_site in site_comparisons.items():
            if test_site['coverage_difference'] > coverage_difference_tolerance:
                out_of_tolerance_sites += 1
        self.assertLessEqual(out_of_tolerance_sites, 50)

    def test_beta_proportion(self):
        # set z threshold
        z_threshold = 3
        # count site with z score above threshold
        z_site_count = 0
        for label, test_site in site_comparisons.items():
            if test_site['beta_z_value'] >= z_threshold:
                z_site_count += 1
        self.assertLessEqual(z_site_count, 250)

    def test_on_target_read_alignments(self):
        # asses proportion of reads that mapped to simulated region
        on_target_alignments = read_stats['On_Prim_PropPair_NoDup'] + read_stats['On_Prim_Discord_NoDup']
        self.assertGreater(on_target_alignments / read_stats['ObservedAlignments'],  0.97)

    def test_off_target_reaad_alignments(self):
        off_target_alignments = read_stats['Off_Prim_PropPair_NoDup'] + read_stats['Off_Prim_Discord_NoDup']
        self.assertLess(off_target_alignments / read_stats['ObservedAlignments'], 0.005)

    def test_unmapped_read_alignments(self):
        self.assertLess(read_stats['UnalignedAlignments'] / read_stats['ObservedAlignments'], 0.001)


if __name__ == '__main__':
    unittest.main()
