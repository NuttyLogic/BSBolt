import random
from typing import Dict, Union
from tqdm import tqdm
from BSBolt.Simulate.SetCyotsineMethylation import SetCytosineMethylation
from BSBolt.Simulate.StreamSim import StreamSim
from BSBolt.Utils.UtilityFunctions import reverse_complement


class SimulateMethylatedReads:
    """Tool to simulated methylation bisulfite sequencing reads. The process works in three distinct steps:
        1. Given a reference Cytosine methylation levels are designated by dinucleotide context, ie CG, CT, etc.
        2. Illumina sequencing Reads are simulated using ART (Huang et al. 2012)
        3. Illumina reads are converted to bisulfite sequencing reads with unconverted methylated cytosines
        Keyword Arguments:
            reference_file (str): path to fasta reference file, all contigs should be in same fasta file
            wgsim_path (str): path to ART executable
     """

    def __init__(self, reference_file: str = None, wgsim_path: str = None, sim_output: str = None,
                 sequencing_error: float = 0.005, mutation_rate: float = 0.0010, mutation_indel_fraction: float = 0.15,
                 indel_extension_probability: float = 0.15, random_seed: int = -1,
                 paired_end: bool = False, read_length: int = 100,
                 read_depth: int = 20, undirectional: bool = False, methylation_reference: str = None,
                 cgmap: str = None, ambiguous_base_cutoff: float = 0.05, haplotype_mode: bool = False,
                 pe_fragment_size: int = 400, insert_deviation: int = 25, mean_insert_size: int = 100,
                 collect_ch_sites: bool = True, collect_sim_stats: bool = False, verbose: bool = True,
                 overwrite_db: bool = False):
        self.sim_command = [wgsim_path, '-1', str(read_length), '-2', str(read_length),
                            '-e', str(sequencing_error), '-d', str(pe_fragment_size),
                            '-s', str(insert_deviation),
                            '-r', str(mutation_rate), '-R', str(mutation_indel_fraction),
                            '-X', str(indel_extension_probability), '-S', str(random_seed),
                            '-A', str(ambiguous_base_cutoff), '-I', str(mean_insert_size)]
        if haplotype_mode:
            self.sim_command.append('-h')
        self.sim_command.append(reference_file)
        self.sim_db = SetCytosineMethylation(reference_file=reference_file,
                                             sim_dir=sim_output,
                                             methylation_reference=methylation_reference,
                                             cgmap=cgmap,
                                             collect_ch_sites=collect_ch_sites,
                                             overwrite_db=overwrite_db)
        self.sim_output = sim_output
        self.paired_end = paired_end
        self.output_objects = self.get_output_objects
        self.undirectional = undirectional
        self.read_coverage = (read_length, read_depth)
        self.reference = self.sim_db.reference
        self.collect_sim_stats = collect_sim_stats
        self.tqdm_disabe = False if verbose else True
        self.current_contig = None
        self.contig_profile = None
        self.contig_values = {}
        self.variant_data = {}

    def run_simulation(self):
        """Commands to execute read simulation"""
        genome_length = sum([len(seq) for seq in self.reference.values()])
        coverage_length = self.read_coverage[0] * 2 if self.paired_end else self.read_coverage[0]
        read_number = (genome_length / coverage_length) * self.read_coverage[1]
        print('Simulating Methylated Reads')
        self.sim_command.extend(['-N', str(int(read_number)), self.sim_db.reference_file])
        self.simulate_methylated_reads()
        print('Finished Simulation')

    def simulate_methylated_reads(self):
        for variant_contig, sim_data in tqdm(StreamSim(paired_end=self.paired_end, sim_command=self.sim_command),
                                             desc='Simulating Bisulfite Converted Read Pairs',
                                             disable=self.tqdm_disabe):
            if variant_contig:
                self.output_reference()
                self.current_contig = variant_contig
                self.contig_profile = self.get_methylation_reference(self.current_contig, sim_data)
                self.variant_data = sim_data
            else:
                if sim_data[1]['chrom'] != self.current_contig:
                    self.output_reference()
                    self.current_contig = sim_data[1]['chrom']
                    self.contig_profile = self.get_methylation_reference(self.current_contig)
                self.process_read_group(sim_data)
        self.sim_db.sim_db.output_contig(self.contig_profile, self.current_contig)
        self.output_reference()
        for output in self.output_objects:
            output.close()

    def output_reference(self):
        self.sim_db.sim_db.output_contig(self.contig_profile, self.current_contig)
        if self.variant_data:
            self.sim_db.sim_db.output_contig(self.variant_data, self.current_contig, variant=True)
        if self.collect_sim_stats:
            self.sim_db.sim_db.output_contig(self.contig_values, self.current_contig, values=True)
            self.contig_values = {}

    def process_read_group(self, sim_data):
        # randomly select reference strand
        sub_pattern = ('C', 'T') if self.random_roll(0.5) else ('G', 'A')
        # set read methylation
        self.set_read_methylation(sim_data[1], sub_base=sub_pattern[0])
        self.set_read_methylation(sim_data[2], sub_base=sub_pattern[0])
        # in silico bisulfite conversion
        sim_data[1]['seq'] = sim_data[1]['seq'].replace(sub_pattern[0], sub_pattern[1]).upper()
        if self.paired_end:
            sim_data[2]['seq'] = sim_data[2]['seq'].replace(sub_pattern[0], sub_pattern[1]).upper()
        # switch subpattern randomly for output if undirectional
        if self.undirectional:
            sub_pattern = ('C', 'T') if self.random_roll(0.5) else ('G', 'A')
        self.output_sim_reads(sim_data, sub_pattern[0])

    def output_sim_reads(self, sim_data, sub_base):
        # format reads
        reverse_read = 2
        if sub_base == 'G':
            reverse_read = 1
        ref_strand = 'W' if sim_data[1]['sub_base'] == 'C' else 'C'
        conversion = 'C2T' if sim_data[1]['sub_base'] == sub_base else 'G2A'
        sim_data[reverse_read]['seq'] = reverse_complement(sim_data[reverse_read]['seq'])
        sim_data[reverse_read]['qual'] = sim_data[reverse_read]['qual'][::-1]
        read_label = f'@{sim_data[1]["read_id"]}/1'
        read_comment = f'+{sim_data[1]["chrom"]}:{sim_data[1]["start"]}:' \
                       f'{sim_data[1]["end"]}:{sim_data[1]["cigar"]}:{ref_strand}{conversion}'
        read = f'{read_label}\n{sim_data[1]["seq"]}\n{read_comment}\n{sim_data[1]["qual"]}\n'
        self.output_objects[0].write(read)
        if self.paired_end:
            read_label = f'@{sim_data[2]["read_id"]}/2'
            read_comment = f'+{sim_data[1]["chrom"]}:{sim_data[2]["start"]}:' \
                           f'{sim_data[2]["end"]}:{sim_data[2]["cigar"]}:{ref_strand}{conversion}'
            read = f'{read_label}\n{sim_data[2]["seq"]}\n{read_comment}\n{sim_data[2]["qual"]}\n'
            self.output_objects[1].write(read)

    def set_read_methylation(self, read, sub_base='C'):
        """Set methylation according to sim value, variants can be methylation but not sequencing errors"""
        ref_seq, seq, cigar = self.reference[read['chrom']][read['start']: read['end']], read['seq'], read['cigar']
        # start iterators
        assert len(seq) == len(cigar)
        start, end = read['start'], read['end']
        methyl_read, methyl_cigar = list(seq), list(cigar)
        methyl_base_info = read['c_base_info'] if sub_base == 'C' else read['g_base_info']
        for site in methyl_base_info.split(','):
            if not site:
                break
            seq_pos, offset = site.split('_')
            seq_pos, offset = int(seq_pos), int(offset)
            cigar_base = cigar[seq_pos]
            ref_pos = start + seq_pos + offset
            seq_base = seq[seq_pos]
            ref_base = ref_seq[seq_pos + offset]
            if cigar_base == 'M':
                # check for sequencing error, if error don't evaluate methylation status
                meth_base, meth_base_cigar = self.handle_match(seq_base, ref_base, ref_pos)
            elif cigar_base == 'X':
                meth_pos = f'{ref_pos}_{ref_base}_{seq_base}'
                meth_base, meth_base_cigar = self.set_variant_methylation(meth_pos, seq_base)
            else:
                meth_pos = f'{ref_pos}_{cigar_base}'
                meth_base, meth_base_cigar = self.set_variant_methylation(meth_pos, seq_base, 'I')
            if meth_base:
                methyl_read[seq_pos] = meth_base
                methyl_cigar[seq_pos] = meth_base_cigar
        read['seq'] = ''.join(methyl_read)
        read['cigar'] = ''.join(methyl_cigar)
        read['sub_base'] = sub_base

    def set_variant_methylation(self, meth_pos, seq_base, variant_type='X'):
        cigar_type = 'Z' if variant_type == 'X' else 'R'
        methyl_status, context = self.set_base_methylation(meth_pos)
        if methyl_status:
            return seq_base.lower(), cigar_type
        else:
            return seq_base, cigar_type.lower()

    def handle_match(self, seq_base, ref_base, pos):
        if seq_base != ref_base:
            return seq_base, 'E'
        else:
            methyl_status, context = self.set_base_methylation(pos)
            if methyl_status:
                return seq_base.lower(), 'C' if context else 'Y'
            else:
                return seq_base, 'c' if context else 'y'

    def set_base_methylation(self, methyl_position):
        # nucleotide, methylation_level, context, 0, 0, 0
        try:
            methyl_info = self.contig_profile[str(methyl_position)]
        except KeyError:
            return False, 0
        methyl_status = self.random_roll(methyl_info[1])
        if self.collect_sim_stats:
            if f'{methyl_position}' not in self.contig_values:
                self.contig_values[f'{methyl_position}'] = [0, 0]
            if methyl_status:
                self.contig_values[f'{methyl_position}'][0] += 1
            else:
                self.contig_values[f'{methyl_position}'][1] += 1
        return methyl_status, methyl_info[2]

    def get_methylation_reference(self, contig: str, variant_data: Union[bool, Dict] = False):
        if variant_data:
            contig_profile = self.sim_db.get_contig_methylation(contig)
            self.sim_db.set_variant_methylation(variant_data, contig_profile, self.current_contig)
            return contig_profile
        else:
            return self.sim_db.get_contig_methylation(contig)

    @property
    def get_output_objects(self):
        output_list = [open(f'{self.sim_output}_1.fq', 'w')]
        if self.paired_end:
            output_list.append(open(f'{self.sim_output}_2.fq', 'w'))
        return output_list

    @staticmethod
    def random_roll(proportion_positive: float = 0.1) -> bool:
        """Return true is random < proportion_positive, used to set individual read methylation"""
        return random.random() < proportion_positive
