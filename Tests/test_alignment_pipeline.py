import gzip
import os
import pickle
import subprocess
import unittest
import numpy as np

# get current directory

test_directory = os.path.dirname(os.path.realpath(__file__))
bsb_directory = '/'.join(test_directory.split('/')[:-1]) + '/'
bsbolt = f'{bsb_directory}BSBolt.py'
# generate simulated reads
bsb_simulate_commands = ['python3', bsbolt, 'Simulate',
                         '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                         '-O', f'{bsb_directory}Tests/TestSimulations/BSB_pe', '-U', '-PE']
subprocess.run(bsb_simulate_commands)

print('Reads Simulated')
# map simulated reads

print('Building Methylation Index')
bsb_index_commands = ['python3', bsbolt, 'Index', '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                      '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB']
subprocess.run(bsb_index_commands)
print('BSB Index Built')


bsb_align_commands = ['python3', bsbolt, 'Align',
                      '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB', '-F1',
                      f'{bsb_directory}Tests/TestSimulations/BSB_pe_meth_1.fastq', '-F2',
                      f'{bsb_directory}Tests/TestSimulations/BSB_pe_meth_2.fastq', '-O',
                      f'{bsb_directory}Tests/BSB_pe_test', '-S', '-OU', '-BT2-k', '10']
subprocess.run(bsb_align_commands)


print('Calling Methylation')

bs_call_methylation_args = ['python3', bsbolt, 'CallMethylation', '-I',
                            f'{bsb_directory}Tests/BSB_pe_test.sorted.bam',
                            '-O', f'{bsb_directory}Tests/BSB_pe_test',
                            '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB',
                            '-t', '6', '-verbose', '-min-qual', '10']
subprocess.run(bs_call_methylation_args)
print('Methylation Values Called')


# import methylation calling dict
all_methylation_sites = {}

output_chromosome = ['chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']
for chrom in output_chromosome:
    with open(f'{bsb_directory}Tests/TestSimulations/BSB_pe.{chrom}.pkl', 'rb') as sim_sites:
        chrom_sites = pickle.load(sim_sites)
        all_methylation_sites.update(chrom_sites['Watson'])
        all_methylation_sites.update(chrom_sites['Crick'])

# import CGmap Calls
cgmap_sites = {}

for line in gzip.open(f'{bsb_directory}Tests/BSB_pe_test.CGmap.gz', 'rb'):
    processed_line = line.decode('utf-8').replace('\n', '').split('\t')
    cgmap_sites[f'{processed_line[0]}:{int(processed_line[2]) - 1}'] = dict(nucleotide=processed_line[1],
                                                                            methylation_level=processed_line[5],
                                                                            context=processed_line[4],
                                                                            methylated_reads=processed_line[6],
                                                                            total_reads=processed_line[7])


def z_test_of_proportion(a_yes, a_no, b_yes, b_no):
    a_total = a_yes + a_no
    b_total = b_yes + b_no
    a_prop = a_yes / a_total
    b_prop = b_yes / b_total
    p_hat = (a_yes + b_yes) / (a_total + b_total)
    try:
        return (a_prop - b_prop) / np.sqrt(p_hat * (1 - p_hat) * (1 / a_total + 1 / b_total))
    except RuntimeWarning:
        return 0

# get line comparisons


site_comparisons = {}
for site, cgmap_values in cgmap_sites.items():
    site_comparison = dict(coverage_difference=0, simulation_beta=0, mapped_beta=0, beta_z_value=0)
    cgmap_site_coverage = int(cgmap_values['total_reads'])
    reference_values = all_methylation_sites[site]
    reference_coverage = int(reference_values['methylated_reads']) + int(reference_values['unmethylated_reads'])
    site_comparison['coverage_difference'] = abs(cgmap_site_coverage - reference_coverage)
    site_comparison['simulation_beta'] = reference_values['methylation_level']
    site_comparison['mapped_beta'] = cgmap_values['methylation_level']
    cgmap_meth = int(cgmap_sites[site]['methylated_reads'])
    cgmap_unmeth = int(cgmap_sites[site]['total_reads']) - int(cgmap_sites[site]['methylated_reads'])
    ref_meth = int(reference_values['methylated_reads'])
    ref_unmeth = int(reference_values['unmethylated_reads'])
    z = abs(z_test_of_proportion(a_yes=cgmap_meth, a_no=cgmap_unmeth, b_yes=ref_meth, b_no=ref_unmeth))
    site_comparison['beta_z_value'] = z
    site_comparisons[site] = site_comparison


class TestBSBPipeline(unittest.TestCase):

    def setUp(self):
        pass

    def test_read_coverage(self):
        # set coverage difference between simulated and mapped reads to consider site out of tolerance
        coverage_difference_tolerance = 5
        # count number of sites out of tolerance
        out_of_tolerance_sites = 0
        for test_site in site_comparisons.values():
            if test_site['coverage_difference'] > coverage_difference_tolerance:
                out_of_tolerance_sites += 1
        self.assertLessEqual(out_of_tolerance_sites, 250)

    def test_beta_proportion(self):
        # set z threshold
        z_threshold = 3
        # count site with z score above threshold
        z_site_count = 0
        for test_site in site_comparisons.values():
            if test_site['beta_z_value'] >= z_threshold:
                z_site_count += 1
        self.assertLessEqual(z_site_count, 250)


if __name__ == '__main__':
    unittest.main()
