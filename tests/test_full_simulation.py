from collections import defaultdict
import unittest
from BSBolt.Simulate import SimulateMethylatedReads
from BSBolt.Utils.AlignmentEvaluation import get_read_reference_info
from BSBolt.Utils.FastqIterator import OpenFastq
from BSBolt.Utils.UtilityFunctions import get_external_paths
from tests.TestHelpers import bsb_directory


bwa_path, wgsim_path = get_external_paths()


sim_out = f'{bsb_directory}tests/TestSimulations/sim_pe'
test_genome = f'{bsb_directory}tests/TestData/BSB_test.fa'


def check_bases(sequence, methyl_cigar, ref_strand, count):
    conversion_check = {'C', 'Y', 'Z', 'R'}
    unconverted_check = {'c', 'y', 'r'}
    conversion_pattern = ('C', 'T')
    if ref_strand[1] == 'G':
        conversion_pattern = ('G', 'A')
    read_check = True
    for base, cigar in zip(sequence, methyl_cigar):
        if cigar in conversion_check:
            if base != conversion_pattern[0]:
                read_check = False
        elif cigar in unconverted_check:
            if base != conversion_pattern[1]:
                read_check = False
    return read_check

# sequencing errors are introduced following methylation setting so set error rate to zero to check cigar string


test_sim = SimulateMethylatedReads(reference_file=test_genome,
                                   paired_end=True,
                                   sim_output=sim_out,
                                   undirectional=True,
                                   collect_sim_stats=True,
                                   overwrite_db=True,
                                   mutation_rate=0.01,
                                   sequencing_error=0.00,
                                   read_length=100,
                                   read_depth=20,
                                   mutation_indel_fraction=0.2)

test_sim.run_simulation()

test_fastqs = [f'{bsb_directory}tests/TestSimulations/sim_pe_1.fq',
               f'{bsb_directory}tests/TestSimulations/sim_pe_2.fq']

sim_data = get_read_reference_info(test_fastqs)

sim_regions = defaultdict(int)


class TestSimulation(unittest.TestCase):

    def setUp(self):
        pass

    def testReferenceSetting(self):
        count = 0
        for fq1_read, fq2_read in zip(*(OpenFastq(test_fastqs[0]), OpenFastq(test_fastqs[1]))):
            count += 1
            fq1_read_info = sim_data[fq1_read[0].replace('@', '')]
            fq2_read_info = sim_data[fq2_read[0].replace('@', '')]
            self.assertEqual(fq1_read[0].split('/')[0], fq2_read[0].split('/')[0])
            self.assertTrue(check_bases(fq1_read[1], fq1_read_info['cigar'], fq1_read_info['ref_strand'], count))
            self.assertTrue(check_bases(fq2_read[1], fq2_read_info['cigar'], fq2_read_info['ref_strand'], count))
            sim_regions[f'{fq1_read_info["chrom"]}:{round(fq1_read_info["start"], -2)}'] += 1
            sim_regions[f'{fq2_read_info["chrom"]}:{round(fq2_read_info["start"], -2)}'] += 1

    def test_coverage(self):
        low_coverage_regions = 0
        high_coverage_regions = 0
        coverage_regions = 0
        # expect 200 read per 100 bp
        for count in sim_regions.values():
            coverage_regions += 1
            if count < 10:
                low_coverage_regions += 1
            elif count > 30:
                high_coverage_regions += 1
        self.assertLess(low_coverage_regions, 500)
        self.assertLess(high_coverage_regions, 500)


if __name__ == '__main__':
    unittest.main()
