import gzip
import os
import pickle
import subprocess
import unittest
import numpy as np

from bsbolt.Utils.CGmapIterator import OpenCGmap
from tests.TestHelpers import bsb_directory, z_test_of_proportion
'''
# generate simulated reads
bsb_simulate_commands = ['python3', '-m', 'bsbolt', 'Simulate',
                         '-G', f'{bsb_directory}tests/TestData/BSB_test.fa',
                         '-O', f'{bsb_directory}tests/TestSimulations/BSB_pe_var', '-PE',
                         '-MR', '0.01', '-MI', '0.0', '-SE', '0.0', '-verbose', '-overwrite', '-RD', '20']
subprocess.run(bsb_simulate_commands)
print('Reads Simulated')
# map simulated reads

if not os.path.exists(f'{bsb_directory}tests/TestData/BSB_Test_DB/BSB_ref.fa.amb'):
    bsb_index_commands = ['python3', '-m', 'bsbolt', 'Index', '-G', f'{bsb_directory}tests/TestData/BSB_test.fa',
                        '-DB', f'{bsb_directory}tests/TestData/BSB_Test_DB']
    subprocess.run(bsb_index_commands)
    print('bsbolt Index Built')


bsb_align_commands = ['python3', '-m', 'bsbolt', 'Align',
                      '-DB', f'{bsb_directory}tests/TestData/BSB_Test_DB', '-F1',
                      f'{bsb_directory}tests/TestSimulations/BSB_pe_var_1.fq', '-F2',
                      f'{bsb_directory}tests/TestSimulations/BSB_pe_var_2.fq', '-O',
                      f'{bsb_directory}tests/BSB_pe_var_test', '-t', '10', '-OT', '2']
subprocess.run(bsb_align_commands)

sorted_output = f'{bsb_directory}tests/BSB_pe_var_test.sorted.bam'
if os.path.exists(sorted_output):
    subprocess.run(['rm', sorted_output])
    subprocess.run(['rm', f'{sorted_output}.bai'])

subprocess.run(['python3', '-m', 'bsbolt', 'Sort', '-I', f'{bsb_directory}tests/BSB_pe_var_test.bam',
                '-O', f'{bsb_directory}tests/BSB_pe_var_test.sorted.bam'])

print('Calling Variation')
'''
subprocess.run(['bsbolt', 'CallVariation', '-I', f'{bsb_directory}tests/BSB_pe_var_test.sorted.bam', '-DB', 
                f'{bsb_directory}tests/TestData/BSB_Test_DB', '-t', '8', '-verbose', '-O',
                f'{bsb_directory}tests/BSB_pe_var_test', '-MQ', '40', '-min', '15'])

# open all chrom output information

variation_db = {}

for contig in ['chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']:
    with open(f'{bsb_directory}tests/TestSimulations/BSB_pe_var.{contig}_variants.pkl', 'rb') as sim:
        sim_var = pickle.load(sim)
        for pos, var in sim_var.items():
            variation_db[f'{var["chrom"]}:{pos}'] = var


var_calls = {}

var_keys = ('contig', 'pos-1', 'pos', 'call_p', 'call_score', 'ref_base', 'genotype')

for line in OpenCGmap(f'{bsb_directory}tests/BSB_pe_var_test.bed.gz'):
    call = {key: value for key, value in zip(var_keys, line)}
    var_calls[f'{call["contig"]}:{call["pos"]}'] = call

false_positives = []
false_negatives = []
true_positives = []
true_negatives = []

called_variants = set()

print(len(var_calls))

for var_id, call in var_calls.items():
    called_genotype = ','.join(sorted(call['genotype'].split(',')))
    try:
        sim_genotype = ','.join(sorted(variation_db[var_id]['iupac']))
    except KeyError:
        sim_genotype = call['ref_base']
    consensus_call = called_genotype == sim_genotype
    if var_id not in variation_db:
        if consensus_call:
            true_positives.append(float[call['call_score']])
        else:
            #print(called_genotype, sim_genotype)
            false_positives.append(float(call['call_score']))
    else:
        if consensus_call:
            true_positives.append(float(call['call_score']))
        else:
            #print(call, variation_db[var_id])
            false_positives.append(float(call['call_score']))

print('end', np.median(np.array(false_positives)), len(false_positives))
print(np.median(np.array(true_positives)), len(true_positives))


class TestBSBPipeline(unittest.TestCase):
    """ The first 5000bp for chr10 are duplicated as chr15 in the simulation reference."""

    def setUp(self):
        pass

    

#if __name__ == '__main__':
#    unittest.main()
