import random
from typing import Dict, List, Union
import numpy as np
from BSBolt.Simulate.SetCyotsineMethylation import SetCytosineMethylation
from BSBolt.Simulate.StreamSim import StreamSim


class SimulateMethylatedReads:
    """Tool to simulated methylation bisulfite sequencing reads. The process works in three distinct steps:
        1. Given a reference Cytosine methylation levels are designated by dinucleotide context, ie CG, CT, etc.
        2. Illumina sequencing Reads are simulated using ART (Huang et al. 2012)
        3. Illumina reads are converted to bisulfite sequencing reads with unconverted methylated cytosines
        Keyword Arguments:
            reference_file (str): path to fasta reference file, all contigs should be in same fasta file
            wgsim_path (str): path to ART executable
            output_path (str): output prefix for processed files
            methylation_reference_output (str): output path for methylation reference, (useful if using built reference)
            paired_end (bool): output paired end reads, default = False / single end reads
            read_length (int): read length
            undirectional (bool): simulate both watson / crick PCR strands
            methylation_reference (str): path to CGmap file to build methylation reference
            methylation_profile (str): previously generated BSBolt simulation database
            insertion_rate1 (float): read 1 insertion rate
            insertion_rate2 (float): read 2 insertion rate
            deletion_rate1 (float): read 1 deletion rate
            deletion_rate2 (float): read 2 deletion rate
            n_base_cutoff (int): threshold for including a read with a gap or ambiguous base
            sequencing_system (str): sequencing system to simulate quality profile
            pe_fragment_size (int): fragment size for simulations
            fragment_size_deviation (int): standard deviation of fragment size for simulation
            read1_quality_profile (str): path to art quality profile for read 1
            read2_quality_profile (str): path to art quality profile for read 2

        Attributes:
            self.simulate_commands (dict): are simulation commands
            self.simulation_out (BSBolt.BSB_Simulation.SimulationOutput): collection of output paths
            self.output_objects (list): list of fastq outputs
            self.paired_end (bool): reads paired
            self.undirectional (bool): simulated pcr products of Watson and Crick strands
            self.reference_contig_size (dict): dict of contig lengths
            self.current_contig (str): label of contig being processed
            self.cpg_distribution (np.vector): binomial distribution to draw CpG methylation proportions
            self.ch_distribution (np.vector): binomial distribution to dra CH methylation proportions
            """

    def __init__(self, reference_file: str = None, wgsim_path: str = None, sim_output: str = None,
                 sequencing_error: float = 0.020, mutation_rate: float = 0.0010, mutation_indel_fraction: float = 0.15,
                 indel_extension_probability: float = 0.15, random_seed: int = -1,
                 methylation_reference_output: str = None, paired_end: bool = False, read_length: int = 125,
                 read_depth: int = 20, undirectional: bool = False, methylation_reference: str = None,
                 methylation_profile: str = None, ambiguous_base_cutoff: float = 0.05, haplotype_mode: bool = False,
                 pe_fragment_size: int = 400, fragment_size_deviation: int = 50,
                 collect_ch_sites: bool = True):
        self.sim_command = [wgsim_path, '-1', str(read_length), '-2', str(read_length),
                            '-e', str(sequencing_error), '-d', str(pe_fragment_size),
                            '-s', str(fragment_size_deviation),
                            '-r', str(mutation_rate), '-R', str(mutation_indel_fraction),
                            '-X', str(indel_extension_probability), '-S', str(random_seed),
                            '-A', str(ambiguous_base_cutoff)]
        if haplotype_mode:
            self.sim_command.append('-h')
        self.sim_db = SetCytosineMethylation(reference_file=reference_file,
                                             sim_dir=sim_output,
                                             methylation_reference=methylation_reference,
                                             cgmap=methylation_profile,
                                             collect_ch_sites=collect_ch_sites)
        self.simulation_out = None
        self.output_path = sim_output
        self.paired_end = paired_end
        self.undirectional = undirectional
        self.read_coverage = (read_length, read_depth)
        self.reference = None
        self.current_contig = None
        self.contig_key = None
        self.contig_values = None
        self.variant_data = {}
        self.output_objects = None

    def run_simulation(self):
        """Commands to execute read simulation"""
        print('Setting Cytosine Methylation')
        genome_length = sum([len(seq) for seq in self.reference.values()])
        coverage_length = self.read_coverage[0] * 2 if self.paired_end else self.read_coverage[0]
        read_number = (genome_length / coverage_length) * self.read_coverage[1]
        print('Simulating Methylated Reads')
        self.sim_command.extend(['-N', str(int(read_number)), self.sim_db.reference_file])

        self.simulate_methylated_reads()
        print('Finished Simulation')

    def simulate_methylated_reads(self):
        for variant_contig, sim_data in StreamSim(paired_end=self.paired_end, sim_command=self.sim_command):
            if variant_contig:
                self.sim_db.sim_db.output_contig(self.contig_key, self.contig_values, self.current_contig)
                self.contig_key, self.contig_values = self.get_methylation_reference(variant_contig, sim_data)
                self.current_contig = variant_contig
                self.variant_data = sim_data
            else:
                if sim_data['chrom'] != self.current_contig:
                    self.sim_db.sim_db.output_contig(self.contig_key, self.contig_values, self.current_contig)
                    self.current_contig = sim_data['chrom']
                    self.contig_key, self.contig_values = self.get_methylation_reference(self.current_contig)
                self.process_read_group(sim_data)

    def process_read_group(self, sim_data):
        undirectional = False
        if self.undirectional:
            undirectional = self.random_roll(0.5)
        subsitution_pattern = 1 if self.random_roll(0.5) else 0

    def set_read_methylation(self, read, sub_nuc='C'):
        ref_seq = iter(self.reference[read['chrom']][read['start']: read['end'] + 1])
        ref_base = next(ref_seq)

    def get_methylation_reference(self, contig: str, variant_data: Union[bool, Dict] = False):
        if variant_data:
            contig_key, contig_values = self.sim_db.get_contig_methylation(contig)
            return self.sim_db.set_variant_methylation(variant_data, contig_key, contig_values)
        else:
            return self.sim_db.get_contig_methylation(contig)

    def get_contig_methylation_reference(self, contig_id):
        if not self.current_contig:
            self.cytosine_dict = self.simulation_out.load_contig(contig_id)
            self.current_contig = contig_id
        elif contig_id != self.current_contig:
            self.simulation_out.output_contig_methylation_reference(self.cytosine_dict, self.current_contig)
            if contig_id != 'stop':
                self.cytosine_dict = self.simulation_out.load_contig(contig_id)
                self.current_contig = contig_id

    @staticmethod
    def nucleotide_check(nucleotide: str) -> bool:
        """Check nucleotide is Cytosine or Guanine"""
        methylation_nucleotides = {'C', 'G'}
        if nucleotide in methylation_nucleotides:
            return True
        return False

    @staticmethod
    def write_fastq(output_object, fastq_line):
        for line in fastq_line:
            output_object.write(f'{line}\n')

    @property
    def get_output_objects(self):
        output_list = [open(f'{self.output_path}_1.fastq', 'w')]
        if self.paired_end:
            output_list.append(open(f'{self.output_path}_2.fastq', 'w'))
        return output_list

    @staticmethod
    def convert_simulated_reads(fastq_lines, methylation_strand, mapping_strand):
        """Convert read sequence based on assigned methylation/ mapping strand"""
        replacement_targets = (('C', 'T'), ('G', 'A'))
        target = 0
        if f'{methylation_strand}_{mapping_strand}' in {'Watson_-', 'Crick_+'}:
            target = 1
        second_line = False
        for fastq_line in fastq_lines:
            # second line preserves methylation strand of first
            if second_line:
                # target will be opposite of first strand
                target = 0 if target == 1 else 1
                fastq_line[1] = fastq_line[1].replace(replacement_targets[target][0],
                                                      replacement_targets[target][1]).upper()
            else:
                fastq_line[1] = fastq_line[1].replace(replacement_targets[target][0],
                                                      replacement_targets[target][1]).upper()
            second_line = True
        return fastq_lines

    @staticmethod
    def random_roll(proportion_positive: float = 0.1) -> bool:
        """Return true is random < proportion_positive, used to set individual read methylation"""
        return random.random() < proportion_positive
