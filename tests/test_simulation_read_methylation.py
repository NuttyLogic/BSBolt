import copy
import unittest
from BSBolt.Simulate.SimulateMethylatedReads import SimulateMethylatedReads
from BSBolt.Utils.UtilityFunctions import reverse_complement, get_external_paths
from tests.TestHelpers import test_directory


bwa_path, wgsim_path, _ = get_external_paths()
# hold read simulation data to test functions


class TestSimOut:

    def __init__(self):
        pass

    @staticmethod
    def write(read):
        return read


sim_out = f'{test_directory}/TestSimulations/wgbs_pe'
test_genome = f'{test_directory}/TestData/BSB_test.fa'

# set rest reference and reads with repeated sequence
test_reference = {'chr10': 'ATCGCATTAA' * 40}
test_read = {1: {'chrom': 'chr10', 'start': 0, 'end': 100,
                 'cigar': 'M' * 100,
                 'pair': 1,
                 'read_id': '1',
                 'comment': '+',
                 'c_base_info': ','.join([f'{2 + 10 * x}_0,{4 + 10 * x}_0' for x in range(10)]),
                 'g_base_info': ','.join([f'{3 + 10 * x}_0' for x in range(10)]),
                 'seq': 'ATCGCATTAA' * 10,
                 'qual': '1' + '2' * 99},
             2: {'chrom': 'chr10', 'start': 100, 'end': 200, 'read_id': '1',
                 'cigar': 'M' * 100,
                 'pair': 2,
                 'comment': '+',
                 'c_base_info': ','.join([f'{2 + 10 * x}_0,{4 + 10 * x}_0' for x in range(10)]),
                 'g_base_info': ','.join([f'{3 + 10 * x}_0' for x in range(10)]),
                 'seq': 'ATCGCATTAA' * 10,
                 'qual': '1' + '2' * 99}
             }
# initialize profile with variant methylation profiles
contig_profile = {'5_A_C': (0, 1, 0),
                  '10_1': (0, 1, 0),
                  '10_2': (1, 1, 0)}
expected_watson_cigar = ['M' for _ in range(400)]
expected_crick_cigar = ['M' for _ in range(400)]

expected_watson_sequence = []
expected_crick_sequence = []

odd = False
for chunk in range(40):
    # chunk CG methylated, CH unmethylated, odd reversed
    if not odd:
        contig_profile[f'{chunk * 10 + 2}'] = (0, 1, 1)
        expected_watson_cigar[chunk * 10 + 2] = 'C'
        contig_profile[f'{chunk * 10 + 3}'] = (1, 1, 1)
        expected_crick_cigar[chunk * 10 + 3] = 'C'
        contig_profile[f'{chunk * 10 + 4}'] = (0, 0, 0)
        expected_watson_cigar[chunk * 10 + 4] = 'y'
        expected_watson_sequence.append('ATCGTATTAA')
        expected_crick_sequence.append('ATCGCATTAA')
        odd = True
    else:
        contig_profile[f'{chunk * 10 + 2}'] = (0, 0, 1)
        expected_watson_cigar[chunk * 10 + 2] = 'c'
        contig_profile[f'{chunk * 10 + 3}'] = (1, 0, 1)
        expected_crick_cigar[chunk * 10 + 3] = 'c'
        contig_profile[f'{chunk * 10 + 4}'] = (0, 1, 0)
        expected_watson_cigar[chunk * 10 + 4] = 'Y'
        expected_watson_sequence.append('ATTGCATTAA')
        expected_crick_sequence.append('ATCACATTAA')
        odd = False


expected_watson_sequence = ''.join(expected_watson_sequence)
expected_crick_sequence = ''.join(expected_crick_sequence)


# all matches
test_read_group_m = copy.deepcopy(test_read)


# mismatch
test_read_group_x = copy.deepcopy(test_read)
seq = list(test_read_group_x[1]['seq'])
seq[5] = 'C'
test_read_group_x[1]['seq'] = ''.join(seq)
cigar = list(test_read_group_x[1]['cigar'])
cigar[5] = 'X'
test_read_group_x[1]['cigar'] = ''.join(cigar)
c_info = test_read_group_x[1]['c_base_info'].split(',')
c_info = c_info[0:2] + ['5_0'] + c_info[2:]
test_read_group_x[1]['c_base_info'] = ','.join(c_info)
x_expected_cigar = list(expected_watson_cigar[0:100])
x_expected_cigar[5] = 'Z'
x_expected_cigar = ''.join(x_expected_cigar)
x_expected_seq = list(expected_watson_sequence)
x_expected_seq[5] = 'C'
x_expected_seq = ''.join(x_expected_seq[0:100])

# deletion
test_read_group_d = copy.deepcopy(test_read)
seq = list(test_read_group_d[1]['seq'])
test_read_group_d[1]['seq'] = ''.join(seq[:10] + seq[12:])
cigar = list(test_read_group_d[1]['cigar'])
test_read_group_d[1]['cigar'] = ''.join(cigar[:10] + cigar[12:])
c_info = test_read_group_d[1]['c_base_info'].split(',')
c_info = c_info[0:2] + [f'{10 * (x + 1)}_2,{2 + 10 * (x +1)}_2' for x in range(9)]
test_read_group_d[1]['c_base_info'] = ','.join(c_info)
d_expected_cigar = list(copy.deepcopy(expected_watson_cigar[0:100]))
d_expected_cigar = ''.join(d_expected_cigar[:10] + d_expected_cigar[12:100])
d_expected_seq = list(expected_watson_sequence)
d_expected_seq = ''.join(d_expected_seq[0:10] + d_expected_seq[12:100])

# insertion
test_read_group_i = copy.deepcopy(test_read)
seq = list(test_read_group_i[1]['seq'])
test_read_group_i[1]['seq'] = ''.join(seq[:10] + ['CG'] + seq[10:])
c_info = test_read_group_i[1]['c_base_info'].split(',')
c_info = c_info[0:2] + ['10_0'] + [f'{4 + 10 * (x + 1)}_-2,{6 + 10 * (x +1)}_-2' for x in range(9)]
test_read_group_i[1]['c_base_info'] = ','.join(c_info)
cigar = list(test_read_group_i[1]['cigar'])
test_read_group_i[1]['cigar'] = ''.join(cigar[:10] + ['12'] + cigar[10:])
i_expected_cigar = list(expected_watson_cigar[0:10] + ['R2'] + expected_watson_cigar[10:100])
i_expected_cigar = ''.join(i_expected_cigar)
i_expected_seq = list(expected_watson_sequence)
i_expected_seq = ''.join(i_expected_seq[0:10] + ['CG'] + i_expected_seq[10:100])

test_sim = SimulateMethylatedReads(reference_file=test_genome,
                                   paired_end=True,
                                   sim_output=sim_out,
                                   undirectional=False,
                                   collect_sim_stats=False)

test_sim.current_contig = 'chr10'
test_sim.reference = test_reference
test_sim.contig_profile = contig_profile
test_sim.output_objects = [TestSimOut(), TestSimOut()]


def output_on_strand(read_group, sub_base='C'):
    while True:
        read_copy = copy.deepcopy(read_group)
        read_output = test_sim.process_read_group(read_copy)
        if read_copy[1]['sub_base'] == sub_base:
            return read_copy, read_output


read_m, output_m = output_on_strand(test_read_group_m, sub_base='C')
read_mc, output_mc = output_on_strand(test_read_group_m, sub_base='G')
read_x, output_x = output_on_strand(test_read_group_x, sub_base='C')
read_d, output_d = output_on_strand(test_read_group_d, sub_base='C')
read_i, output_i = output_on_strand(test_read_group_i, sub_base='C')


class TestReadMethylationSetting(unittest.TestCase):

    def setUp(self):
        pass

    def test_read_m_w_seq(self):
        self.assertEqual(read_m[1]['seq'], expected_watson_sequence[0:100])

    def test_read_m_w_cigar(self):
        self.assertEqual(read_m[1]['cigar'], ''.join(expected_watson_cigar[0:100]))

    def test_read_m_w_2_seq(self):
        self.assertEqual(reverse_complement(read_m[2]['seq']), ''.join(expected_watson_sequence[100:200]))

    def test_read_m_w_2_cigar(self):
        self.assertEqual(read_m[2]['cigar'], ''.join(expected_watson_cigar[100:200:])[::-1])

    def test_crick_seq(self):
        self.assertEqual(read_mc[1]['seq'], reverse_complement(expected_crick_sequence[0:100]))

    def test_crick_cigar(self):
        self.assertEqual(read_mc[1]['cigar'], ''.join(expected_crick_cigar[0:100])[::-1])

    def test_crick_seq_2(self):
        self.assertEqual(read_mc[2]['seq'], ''.join(expected_crick_sequence[100:200]))

    def test_crick_cigar_2(self):
        self.assertEqual(read_mc[2]['cigar'], ''.join(expected_crick_cigar[100:200]))

    def test_read_x_seq(self):
        self.assertEqual(read_x[1]['seq'], x_expected_seq)

    def test_read_x_cigar(self):
        self.assertEqual(read_x[1]['cigar'], x_expected_cigar)

    def test_read_d_seq(self):
        self.assertEqual(read_d[1]['seq'], d_expected_seq)

    def test_read_d_cigar(self):
        self.assertEqual(read_d[1]['cigar'], d_expected_cigar)

    def test_read_i_seq(self):
        self.assertEqual(read_i[1]['seq'], i_expected_seq)

    def test_read_i_cigar(self):
        self.assertEqual(read_i[1]['cigar'], i_expected_cigar)


if __name__ == '__main__':
    unittest.main()


