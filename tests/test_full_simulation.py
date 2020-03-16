import unittest
from BSBolt.Simulate.SimulateMethylatedReads import SimulateMethylatedReads
from BSBolt.Utils.FastqIterator import OpenFastq
from BSBolt.Utils.UtilityFunctions import get_external_paths
from tests.TestHelpers import test_directory


bwa_path, wgsim_path = get_external_paths()


sim_out = f'{test_directory}/TestSimulations/wgbs_pe'
test_genome = f'{test_directory}/TestData/BSB_test.fa'


test_sim = SimulateMethylatedReads(reference_file=test_genome,
                                   wgsim_path=wgsim_path,
                                   paired_end=True,
                                   sim_output=sim_out,
                                   undirectional=False,
                                   collect_sim_stats=False,
                                   overwrite_db=True,
                                   mutation_rate=0.0,
                                   sequencing_error=0.0,
                                   read_length=100,
                                   read_depth=20)

test_sim.simulate_methylated_reads()

# import reads

for read in OpenFastq(f'{sim_out}_1.fq'):
    pass #print(read)



class TestSimulation(unittest.TestCase):

    def setUp(self):
        pass




if __name__ == '__main__':
    unittest.main()
