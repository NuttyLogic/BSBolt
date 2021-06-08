import os
import subprocess
import unittest
import pysam

from bsbolt.Utils.AlignmentEvaluation import get_read_reference_info
from tests.TestHelpers import bsb_directory


bsb_wgbs_masked_index_commands = ['python3', '-m', 'bsbolt', 'Index', '-G',
                                  f'{bsb_directory}tests/TestData/BSB_test.fa',
                                  '-DB', f'{bsb_directory}tests/TestData/BSB_Test_DB_wgbs_masked', '-MR',
                                  f'{bsb_directory}tests/TestData/test_wgbs_masking.bed']

subprocess.run(bsb_wgbs_masked_index_commands)


bsb_simulate_commands = ['python3', '-m', 'bsbolt', 'Simulate', '-CH', '-NF', '0.001',
                         '-G', f'{bsb_directory}tests/TestData/BSB_Test_DB_wgbs_masked/BSB_ref.fa',
                         '-O', f'{bsb_directory}tests/TestSimulations/BSB_pe_masked', '-PE', '-RD', '1',
                         '-verbose']
subprocess.run(bsb_simulate_commands)


bsb_align_commands = ['python3', '-m', 'bsbolt', 'Align',
                      '-DB', f'{bsb_directory}tests/TestData/BSB_Test_DB_wgbs_masked', '-F1',
                      f'{bsb_directory}tests/TestSimulations/BSB_pe_masked_1.fq', '-F2',
                      f'{bsb_directory}tests/TestSimulations/BSB_pe_masked_2.fq', '-O',
                      f'{bsb_directory}tests/BSB_pe_test_masked']
subprocess.run(bsb_align_commands)

sorted_output = f'{bsb_directory}tests/BSB_pe_test_masked.sorted.bam'
if os.path.exists(sorted_output):
    subprocess.run(['rm', sorted_output])

subprocess.run(['python3', '-m', 'bsbolt', 'Sort', '-I', f'{bsb_directory}tests/BSB_pe_test_masked.bam',
                '-O', f'{bsb_directory}tests/BSB_pe_test_masked.sorted.bam'])

if os.path.exists(f'{sorted_output}.bai'):
    subprocess.run(['rm', f'{sorted_output}.bai'])

subprocess.run(['python3', '-m', 'bsbolt', 'BamIndex', '-I', sorted_output])

# import alignment regions

target_regions = []

with open(f'{bsb_directory}tests/TestData/test_wgbs_masking.bed', 'r') as bed:
    for line in bed:
        chrom, start, end = line.replace('\n', '').split('\t')
        start, end = int(start), int(end)
        target_regions.append([chrom, start, end])

target_regions.sort(key=lambda x: x[1])
target_regions.sort(key=lambda x: x[0])

combined_target_regions = {chrom: [] for chrom in set([x[0] for x in target_regions])}

region_chrom, region_start, region_end = target_regions[0][0], target_regions[0][1], target_regions[0][2]
for region in target_regions[1:]:
    if region[0] != region_chrom:
        combined_target_regions[region_chrom].append([region_start, region_end])
        region_chrom, region_start, region_end = region
    else:
        if region_start <= region[1] <= region_end:
            region_end = region[2]
        else:
            combined_target_regions[region_chrom].append([region_start, region_end])
            region_chrom, region_start, region_end = region
combined_target_regions[region_chrom].append([region_start, region_end])

# retrieve read names for reads that overlap target regions


def check_overlap(region, chrom, start, end):
    if chrom != region[0]:
        return False
    if region[1] <= start and end <= region[2]:
        return True
    return False


reference_info = get_read_reference_info([f'{bsb_directory}tests/TestSimulations/BSB_pe_masked_1.fq',
                                          f'{bsb_directory}tests/TestSimulations/BSB_pe_masked_2.fq'])

overlapping_reads = []

for name, info in reference_info.items():
    for site in target_regions:
        if check_overlap(site, info['chrom'], int(info['start']), int(info['end'])):
            overlapping_reads.append(name)

# extract mapped read names

on_target = []
off_target = []

masked_alingnment_file = pysam.AlignmentFile(f'{bsb_directory}tests/BSB_pe_test_masked.sorted.bam', 'rb')
for count, read in enumerate(masked_alingnment_file.fetch()):
    read_name = f'{read.query_name}/1'
    if not read.is_read2:
        read_name = f'{read.query_name}/2'
    if read_name in overlapping_reads:
        on_target.append(read_name)
    else:
        off_target.append(read_name)


class TestMaskedAlignment(unittest.TestCase):

    def setUp(self):
        pass

    def test_reads_on_target(self):
        on_target_proportion = len(on_target) / len(overlapping_reads)
        self.assertGreater(on_target_proportion, .70)

    def test_off_target_reads(self):
        off_target_proportion = len(off_target) / len(overlapping_reads)
        self.assertLess(off_target_proportion, 2)


if __name__ == '__main__':
    unittest.main()
