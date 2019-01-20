import gzip
import os
import pickle
import subprocess
import unittest
import numpy as np

# get current directory

test_directory = os.path.dirname(os.path.realpath(__file__))
bsb_directory = '/'.join(test_directory.split('/')[:-1]) + '/'
bowtie2_path = 'bowtie2'
art_path = '/Users/colinfarrell/Downloads/art_bin_MountRainier/art_illumina'
print('Generating Simulated Methylation Reads')
# generate simulated reads
bsb_simulate_commands = ['python3', f'{bsb_directory}BSBolt-Simulate.py',
                         '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                         '-A', art_path, '-O', f'{bsb_directory}Tests/TestSimulations/BSB_pe', '-U', '-PE']
subprocess.run(bsb_simulate_commands)

print('Reads Simulated')
# map simulated reads

print('Building Methylation Index')
bsb_index_commands = ['python3', f'{bsb_directory}BSBolt-Index.py', '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                      '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB', '-BT2', bowtie2_path]
subprocess.run(bsb_index_commands)
print('BSB Index Built')


bsb_align_commands = ['python3', f'{bsb_directory}BSBolt-Align.py',
                      '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB', '-BT2', bowtie2_path, '-F1',
                      f'{bsb_directory}Tests/TestSimulations/BSB_pe_meth_1.fastq', '-F2',
                      f'{bsb_directory}Tests/TestSimulations/BSB_pe_meth_2.fastq', '-O',
                      f'{bsb_directory}Tests/BSB_pe_test', '-S']
subprocess.run(bsb_align_commands)


print('Calling Methylation')

bs_call_methylation_args = ['python3', f'{bsb_directory}BSBolt-Call-Methylation.py', '-I',
                            f'{bsb_directory}Tests/BSB_pe_test.sorted.bam',
                            '-O', f'{bsb_directory}Tests/BSB_pe_test',
                            '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB',
                            '-t', '4', '-verbose']
subprocess.run(bs_call_methylation_args)
print('Methylation Values Called')


# import methylation calling dict

with open(f'{bsb_directory}Tests/TestSimulations/BSB_pe_methylation_value_dict.pkl', 'rb') \
        as sim_sites:
    methylation_sites = pickle.load(sim_sites)

# import CGmap Calls
cgmap_sites = {}

for line in gzip.open(f'{bsb_directory}Tests/BSB_pe_test.CGmap.gz', 'rb'):
    processed_line = line.decode('utf-8').replace('\n', '').split('\t')
    cgmap_sites[f'{processed_line[0]}:{int(processed_line[2]) - 1}'] = dict(nucleotide=processed_line[1],
                                                                            methylation_level=processed_line[5],
                                                                            context=processed_line[4],
                                                                            methylated_reads=processed_line[6],
                                                                            total_reads=processed_line[7])

all_methylation_sites = dict(methylation_sites['Watson'])
all_methylation_sites.update(methylation_sites['Crick'])
del methylation_sites


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
        for site in site_comparisons.values():
            if site['coverage_difference'] > coverage_difference_tolerance:
                out_of_tolerance_sites += 1
        self.assertLessEqual(out_of_tolerance_sites, 50)

    def test_beta_proportion(self):
        # set z threshold
        z = 3
        # count site with z score above threshold
        z_site_count = 0
        for site in site_comparisons.values():
            if site['beta_z_value'] >= z:
                print(z)
                z_site_count += 1
        self.assertLessEqual(z_site_count, 50)


if __name__ == '__main__':
    unittest.main()
