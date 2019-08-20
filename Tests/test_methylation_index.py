import gzip
import os
import subprocess
import unittest
from BSB.BSB_Index.ProcessCutSites import ProcessCutSites
from BSB.BSB_Utils.UtilityFunctions import reverse_complement

# get current directory

test_directory = os.path.dirname(os.path.realpath(__file__))
bsb_directory = '/'.join(test_directory.split('/')[:-1]) + '/'
bsbolt = f'{bsb_directory}BSBolt.py'

# generate methylation indices for rrbs, wgbs, and wgbs masked

bsb_wgbs_index_commands = ['python3', bsbolt, 'Index', '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                           '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB_wgbs']

subprocess.run(bsb_wgbs_index_commands)

bsb_wgbs_masked_index_commands = ['python3', bsbolt, 'Index', '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                                  '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB_wgbs_masked', '-MR',
                                  f'{bsb_directory}Tests/TestData/test_wgbs_masking.bed']

subprocess.run(bsb_wgbs_masked_index_commands)


bsb_rrbs_index_commands = ['python3', bsbolt, 'Index', '-G', f'{bsb_directory}Tests/TestData/BSB_test.fa',
                           '-DB', f'{bsb_directory}Tests/TestData/BSB_Test_DB_rrbs', '-rrbs', '-MR',
                           f'{bsb_directory}Tests/TestData/test_wgbs_masking.bed', '-rrbs-cut-format', 'CAG-G']

subprocess.run(bsb_rrbs_index_commands)

# get unmasked sequence ranges

mappable_test_regions = []

with open(f'{bsb_directory}Tests/TestData/BSB_Test_DB_wgbs_masked/BSB_ref.fa', 'r') as bsb_ref:
    chrome: str = None
    start: int = None
    end: int = None
    for line in bsb_ref:
        if '>' in line:
            chrome = line.replace('\n', '').replace('>', '')
        else:
            for count, base in enumerate(line.replace('\n', '')):
                if base != '-':
                    if start:
                        end = count
                    else:
                        start = count
                else:
                    if start:
                        mappable_test_regions.append(f'{chrome}:{start}-{end}')
                        start = None
                        end = None

watson_regions = []
crick_regions = []

for region in mappable_test_regions:
    pos = region.split(':')[-1]
    start, end = [int(location) for location in pos.split('-')]
    if 'crick' in region:
        crick_regions.append(abs(start - end))
    else:
        watson_regions.append(abs(start - end))


mappable_bed = []

with open(f'{bsb_directory}Tests/TestData/test_wgbs_masking.bed', 'r') as bed:
    for line in bed:
        line_split = line.replace('\n', '').split('\t')
        mappable_bed.append(f'{line_split[0]}:{line_split[1]}-{line_split[2]}')

mappable_test_regions.sort()
mappable_bed.sort()

# get mappable rrbs regions:
mappable_rrbs_regions = []

with gzip.open(f'{bsb_directory}Tests/TestData/BSB_Test_DB_rrbs/mappable_regions.txt.gz') as mappable_rrbs:
    for line in mappable_rrbs:
        line = line.decode('UTF-8')
        chrom, start, end, seq = line.split('\t')
        mappable_rrbs_regions.append([chrom, int(start), int(end)])


reference_sequences = {}
# retrieve contig sequences
with open(f'{bsb_directory}Tests/TestData/BSB_Test_DB_rrbs/BSB_ref.fa', 'r') as bsb_ref:
    chrome: str = None
    for line in bsb_ref:
        if '>' in line:
            chrome = line.replace('\n', '').replace('>', '')
        else:
            reference_sequences[chrome] = line.replace('\n', '')

reference_sequences_restricted_mapping = {}
# retrieve contig sequences
with open(f'{bsb_directory}Tests/TestData/BSB_Test_DB_wgbs_masked/BSB_ref.fa', 'r') as bsb_ref:
    chrome: str = None
    for line in bsb_ref:
        if '>' in line:
            chrome = line.replace('\n', '').replace('>', '')
        else:
            if 'crick' in chrome:
                reference_sequences_restricted_mapping[chrome] = reverse_complement(line.replace('\n', ''))
            else:
                reference_sequences_restricted_mapping[chrome] = line.replace('\n', '')


cut_site_offsets = list(ProcessCutSites(cut_format='C-CGG,CA-TG,CAG-G').restriction_site_dict.values())
cut_site_sequence = list(ProcessCutSites(cut_format='C-CGG,CA-TG,CAG-G').restriction_site_dict.keys())


class TestReadSimulation(unittest.TestCase):

    def setUp(self):
        pass

    def test_mappable_regions(self):
        # assert mappable regions are actually exposed in reference file
        for map_region in mappable_bed:
            self.assertIn(map_region, mappable_test_regions)

    def test_mappable_size(self):
        # test mappable regions for watson and crick strands are the same size
        for map_region in watson_regions:
            self.assertIn(map_region, crick_regions)

    def test_mappable_sequence(self):
        for chromosome in ['chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']:
            self.assertEqual(reference_sequences_restricted_mapping[chromosome].replace('G', 'A'),
                             reference_sequences_restricted_mapping[f'{chromosome}_crick_bs'].replace('C', 'T'))

    def test_site_offsets(self):
        self.assertEqual(cut_site_offsets, [1, 2, 3, 1])

    def test_cut_sequnece(self):
        self.assertEqual(cut_site_sequence, ['CCGG', 'CATG', 'CAGG', 'CCTG'])

    def test_rrbs_mappable_regions(self):
        for map_region in mappable_rrbs_regions:
            sequence = reference_sequences[map_region[0]][map_region[1]:map_region[2] + 1]
            self.assertTrue('-' not in sequence)

    def test_rrbs_mappable_sites(self):
        for region in mappable_rrbs_regions:
            sequence = reference_sequences[region[0]][region[1]:region[2] + 1]
            seq_start = sequence[0:3]
            seq_end = sequence[-3:]
            if seq_start[0] != 'G':
                self.assertEqual(seq_start, 'TTG')
            if seq_end[-1] != 'T':
                self.assertEqual(seq_end, 'TAG')


if __name__ == '__main__':
    unittest.main()
