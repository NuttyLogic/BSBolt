import gzip
import os
import subprocess
import unittest
import pysam

# get current directory

test_directory = os.path.dirname(os.path.realpath(__file__))
bsb_directory = '/'.join(test_directory.split('/')[:-1]) + '/'
bsbolt = f'{bsb_directory}BSBolt.py'


bsb_simulate_commands = ['python3', bsbolt, 'Simulate',
                         '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                         '-O', f'{bsb_directory}Tests/TestSimulations/BSB_pe', '-U', '-PE']
#subprocess.run(bsb_simulate_commands)

print('Reads Simulated')


bsb_wgbs_masked_index_commands = ['python3', bsbolt, 'Index', '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                                  '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB_wgbs_masked', '-MR',
                                  f'{bsb_directory}Tests/TestData/test_wgbs_masking.bed']

subprocess.run(bsb_wgbs_masked_index_commands)

bsb_align_commands = ['python3', bsbolt, 'Align',
                      '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB_wgbs_masked', '-F1',
                      f'{bsb_directory}Tests/TestSimulations/BSB_pe_meth_1.fastq', '-F2',
                      f'{bsb_directory}Tests/TestSimulations/BSB_pe_meth_2.fastq', '-O',
                      f'{bsb_directory}Tests/BSB_pe_test_masked', '-S']
subprocess.run(bsb_align_commands)

# import alignment regions

target_regions = []

with open(f'{bsb_directory}Tests/TestData/test_wgbs_masking.bed', 'r') as bed:
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


def check_overlap(region, start, end):
    if region[0] <= start <= region[1] and region[0] <= end <= region[1]:
        return True
    return False


overlaping_reads = []

with open(f'{bsb_directory}Tests/TestSimulations/BSB_pe.sam', 'r') as sam:
    sam_iter = iter(sam)
    while True:
        try:
            sam1 = next(sam_iter).split('\t')
        except StopIteration:
            break
        else:
            if '@' in sam1[0]:
                continue
        sam2 = next(sam_iter).split('\t')
        s1_name, s1_chrom, s1_start, s1_end = sam1[0], sam1[2], int(sam1[3]), int(sam1[3]) + 125
        s2_name, s2_chrom, s2_start, s2_end = sam2[0], sam2[2], int(sam2[3]), int(sam2[3]) + 125
        assert s1_name == s2_name
        target_regions = combined_target_regions.get(s1_chrom, [])
        complete_overlap = False
        for site in target_regions:
            if check_overlap(site, s1_start, s1_end) and check_overlap(site, s2_start, s2_end):
                overlaping_reads.append(s1_name)


# extract mapped read names

expected_reads = []
non_target_reads = []

masked_alingment_file = pysam.AlignmentFile(f'{bsb_directory}Tests/BSB_pe_test_masked.sorted.bam', 'rb')
for count, read in enumerate(masked_alingment_file.fetch()):
    if read.query_name in overlaping_reads:
        if read.query_name not in expected_reads:
            expected_reads.append(read.query_name)
    else:
        if read.query_name not in non_target_reads:
            non_target_reads.append(read.query_name)

print(len(expected_reads))
print(len(non_target_reads))
print(len(overlaping_reads))
print(sorted([read for read in overlaping_reads if read not in expected_reads]))
print(sorted(non_target_reads))


