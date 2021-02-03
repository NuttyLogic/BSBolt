import unittest
from bsbolt.Simulate.SetCyotsineMethylation import SetCytosineMethylation
from bsbolt.Utils.UtilityFunctions import get_external_paths
from tests.TestHelpers import test_directory


bwa_path, wgsim_path, _ = get_external_paths()
# hold read simulation data to test functions

sim_out = f'{test_directory}/TestSimulations/wgbs_pe/'
test_genome = f'{test_directory}/TestData/BSB_test.fa'

# set rest reference and reads with repeated sequence
test_reference = {'chr10': 'ATCGCATTAA' * 40}

test_setter = SetCytosineMethylation(reference_file=test_genome,
                                     sim_dir=sim_out,
                                     collect_ch_sites=True,
                                     overwrite_db=True)
test_setter.reference = test_reference
test_profile = test_setter.get_contig_methylation('chr10')


class TestReferenceSetting(unittest.TestCase):

    def setUp(self):
        pass

    def test_cg_c(self):
        for c_cg in range(40):
            site_profile = test_profile[f'{c_cg * 10 + 2}']
            self.assertEqual(site_profile[0], 0)
            self.assertEqual(site_profile[2], 1)

    def test_cg_g(self):
        for g_cg in range(40):
            site_profile = test_profile[f'{g_cg * 10 + 3}']
            self.assertEqual(site_profile[0], 1)
            self.assertEqual(site_profile[2], 1)

    def test_ch(self):
        for ch in range(40):
            site_profile = test_profile[f'{ch * 10 + 4}']
            self.assertEqual(site_profile[0], 0)
            self.assertEqual(site_profile[2], 0)


if __name__ == '__main__':
    unittest.main()
