import os
import unittest
from BSBolt.Simulate.SimulateMethylatedReads import SimulateMethylatedReads
from BSBolt.Utils.UtilityFunctions import reverse_complement, get_external_paths
from tests.TestHelpers import test_directory


bwa_path, wgsim_path = get_external_paths()
# hold read simulation data to test functions


simulation_files = f'{test_directory}/TestSimulations/wgbs_pe'
test_genome = f'{test_directory}/TestData/BSB_test.fa'

test_simulation_files = SimulateMethylatedReads(reference_file=test_genome,
                                                wgsim_path=wgsim_path,
                                                paired_end=True,
                                                sim_output=simulation_files,
                                                undirectional=True,
                                                collect_sim_stats=True)
test_simulation_files.run_simulation()



class TestReadSimulation(unittest.TestCase):

    def setUp(self):
        pass



if __name__ == '__main__':
    unittest.main()

