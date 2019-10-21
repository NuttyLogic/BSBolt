import os
import pickle
import queue
import time
import unittest
import numpy as np
from BSBolt.CallMethylation.CallMethylationValues import CallMethylationValues
from BSBolt.CallMethylation.CallMethylationVector import CallMethylationVector

# set test_directories
# tests will only run if data has previously been simulated

test_directory = os.path.dirname(os.path.realpath(__file__))
bsb_directory = '/'.join(test_directory.split('/')[:-1]) + '/'


bsb_index = os.path.exists(f'{bsb_directory}tests/TestSimulations/BSB_pe.chr11.pkl')
bsb_alignment_file = os.path.exists(f'{bsb_directory}tests/BSB_pe_test.sorted.bam')

if not all([bsb_index, bsb_alignment_file]):
    print('Simulate date before running tests')
    raise ValueError

chr11_methylation_sites = {}

with open(f'{bsb_directory}tests/TestSimulations/BSB_pe.chr11.pkl', 'rb') as sim_sites:
    sim_values = pickle.load(sim_sites)
    chr11_methylation_sites.update(sim_values['Watson'])
    chr11_methylation_sites.update(sim_values['Crick'])


# call methylation values

value_queue = queue.Queue()

value_start = time.time()

value_caller = CallMethylationValues(input_file=f'{bsb_directory}tests/BSB_pe_test.sorted.bam',
                                     genome_database=f'{bsb_directory}tests/TestData/BSB_Test_DB/',
                                     contig='chr11', return_queue=value_queue, ignore_orphans=False,
                                     ignore_overlap=True, min_base_quality=1)
value_caller.call_methylation()

print(f'Values Called: {time.time() - value_start: .3f} sec')

all_values = []


while True:
    values = value_queue.get()
    all_values.extend(values)
    if value_queue.empty():
        break

# process methylation values

cleaned_values = {}

for value in all_values:
    site_label = f'{value[8]}:{value[7] - 1}'
    cleaned_values[site_label] = [value[1], value[2]]


vector_queue = queue.Queue()

vector_start = time.time()

vector_caller = CallMethylationVector(input_file=f'{bsb_directory}tests/BSB_pe_test.sorted.bam',
                                      genome_database=f'{bsb_directory}tests/TestData/BSB_Test_DB/',
                                      contig='chr11', return_queue=vector_queue)
vector_caller.call_methylation()

print(f'Vectors Called: {time.time() - vector_start: .3f} sec')


all_vectors = []


while True:
    values = vector_queue.get()
    all_vectors.extend(values[1])
    if vector_queue.empty():
        break

vector_point_values = {}

vector_point_start = time.time()


def get_vector_point_values(vector, chromosome, value_dict):
    for meth_call, pos in zip(vector[3], vector[4]):
        if f'{chromosome}:{pos}' in value_dict:
            if meth_call:
                value_dict[f'{chromosome}:{pos}'][0] += 1
            else:
                value_dict[f'{chromosome}:{pos}'][1] += 1
        else:
            if meth_call:
                value_dict[f'{chromosome}:{pos}'] = [1, 0]
            else:
                value_dict[f'{chromosome}:{pos}'] = [0, 1]


# get point values
for vector in all_vectors:
    get_vector_point_values(vector, 'chr11', vector_point_values)

print(f'Vector Point Values Called: {time.time() - vector_point_start: .3f} sec')

vector_point_value_comparison = {}

mismatch_count = 0
differences = []
for site, vector_values in vector_point_values.items():
    pileup_values = cleaned_values[site]
    if vector_values != pileup_values:
        mismatch_count += 1
        print(vector_values, pileup_values, chr11_methylation_sites[site], site)
        difference = abs(vector_values[0] - pileup_values[0]) + abs(vector_values[1] - pileup_values[1])
        differences.append(difference)

print(max(differences), np.mean(differences))
print(mismatch_count, len(vector_point_values))


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


for site, site_values in chr11_methylation_sites.items():
    site_comparison = dict(coverage_difference=0, simulation_beta=0, mapped_beta=0, beta_z_value=0)




