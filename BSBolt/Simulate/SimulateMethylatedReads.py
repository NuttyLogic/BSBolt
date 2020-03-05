import random
import subprocess
from typing import Dict, List
from BSBolt.Utils.FastqIterator import OpenFastq
from BSBolt.Utils.AlnIterator import OpenAln
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

    def __init__(self, reference_file: str = None, wgsim_path: str = None, output_path: str = None,
                 sequencing_error: float = 0.020, mutation_rate: float = 0.0010, mutation_indel_fraction: float = 0.15,
                 indel_extension_probability: float = 0.15, random_seed: int = -1,
                 methylation_reference_output: str = None, paired_end: bool = False, read_length: int = 125,
                 read_depth: int = 20, undirectional: bool = False, methylation_reference: str = None,
                 methylation_profile: str = None, ambiguous_base_cutoff: float = 0.05, haplotype_mode: bool = False,
                 pe_fragment_size: int = 400, fragment_size_deviation: int = 50,
                 genome_length: int = None):
        self.simulate_commands = [wgsim_path, '-1', str(read_length), '-2', str(read_length),
                                  '-e', str(sequencing_error), '-d', str(pe_fragment_size),
                                  '-s', str(fragment_size_deviation),
                                  '-r', str(mutation_rate), '-R', str(mutation_indel_fraction),
                                  '-X', str(indel_extension_probability), '-S', str(random_seed),
                                  '-A', str(ambiguous_base_cutoff)]
        if haplotype_mode:
            self.simulate_commands.append('-h')
        if not methylation_reference_output:
            methylation_reference_output = output_path
        self.cytosine_dict_profile_kwargs = dict(reference_file=reference_file,
                                                 methylation_reference_output=methylation_reference_output,
                                                 methylation_reference=methylation_reference,
                                                 methylation_profile=methylation_profile)
        self.simulation_out = None
        self.output_path = output_path
        self.paired_end = paired_end
        self.undirectional = undirectional
        self.read_coverage = (read_length, read_depth)
        self.reference_contig_size = None
        self.current_contig = None
        self.cytosine_dict = None
        self.output_objects = None

    def run_simulation(self):
        """Commands to execute read simulation"""
        print('Setting Cytosine Methylation')
        cytosine_profile = SetCytosineMethylation(**self.cytosine_dict_profile_kwargs)
        self.reference_contig_size, self.simulation_out = cytosine_profile.set_simulated_methylation()
        genome_length = sum([len(seq) for seq in self.reference_contig_size.values()])
        coverage_length = self.read_coverage[0] * 2 if self.paired_end else self.read_coverage[0]
        read_number = (genome_length / coverage_length) * self.read_coverage[1]
        print('Simulating Methylated Reads')
        self.simulate_commands.extend(['-N', str(int(read_number)),
                                       self.cytosine_dict_profile_kwargs['reference_file']])
        read_stream = StreamSim(paired_end=self.paired_end)
        read_stream.stream_sim_reads(self.simulate_commands)
        self.output_objects = self.get_output_objects
        print('Finished Simulation')

    def simulate_illumina_reads(self):
        """Launch external ART command to simulate illumina reads. Insertion and deletion rate set to zero
        to simplify read simulation """
        art = subprocess.run(args=self.simulate_commands)
        if art.returncode:
            print('ART not correctly initialized, please check path')
            raise OSError

    def simulate_methylated_reads(self, aln_files: List[str], fastq_files: List[str]):
        simulation_iterators = [OpenAln(aln_file) for aln_file in aln_files]
        # noinspection PyTypeChecker
        simulation_iterators.extend([OpenFastq(fastq) for fastq in fastq_files])
        for line in zip(*simulation_iterators):
            aln1_profile = self.parse_aln_line(line[0])
            methylation_strand = 'Watson'
            # if undirectional strand assigned regardless of strand else Watson if reference simulation strand positive
            if self.undirectional:
                if self.random_roll(proportion_positive=0.5):
                    methylation_strand = 'Crick'
            else:
                if aln1_profile['reference_strand'] == '-':
                    methylation_strand = 'Crick'
            methylated_reads = []
            if self.paired_end:
                aln2_profile = self.parse_aln_line(line[1])
                methylated_reads.append(self.set_simulated_methylation(aln1_profile, line[2], methylation_strand))
                methylated_reads.append(self.set_simulated_methylation(aln2_profile, line[3], methylation_strand))
            else:
                methylated_reads.append(self.set_simulated_methylation(aln1_profile, line[1], methylation_strand))
            processed_reads = self.convert_simulated_reads(methylated_reads,
                                                           methylation_strand,
                                                           aln1_profile['reference_strand'])
            for output, read in zip(self.output_objects, processed_reads):
                self.write_fastq(output, read)
        self.get_contig_methylation_reference('stop')
        # close fastq output objects
        for output in self.output_objects:
            output.close()

    def set_simulated_methylation(self, aln_profile: Dict, fastq_line: List, methylation_strand: str) -> List:
        """Simulated methylated bases are converted based on case. Uppercase bases are converted lower case bases
         are ignored. Methylation of bases is set based on the cytosine dict.
         Keyword Arguments:
             aln_profile (dict): processed aln file line
             fastq_line (list): list of strings, 4 per fastq line
             methylation_strand (str): strand to use for setting methylation values.
        Return:
            fatst_line (list): fastq line with altered read sequence, (methylated bases set to lower case)
        """
        genome_position = aln_profile['read_index']
        self.get_contig_methylation_reference(aln_profile['read_contig'])
        strand_cytosine: dict = self.cytosine_dict[methylation_strand]
        fastq_read = []
        for reference_nuc, read_nuc in zip(aln_profile['reference_sequence'], aln_profile['read_sequence']):
            # if reference nuc indicates insertion in read skip base and go to next nuc without changing genome_position
            if reference_nuc == '-':
                fastq_read.append(read_nuc)
                continue
            # proceed with processing if both nucleotides either a C or a G
            elif self.nucleotide_check(reference_nuc) and self.nucleotide_check(read_nuc):
                location_meth_info: dict = strand_cytosine.get(f'{aln_profile["read_contig"]}:{genome_position}', False)
                if location_meth_info:
                    if self.random_roll(proportion_positive=location_meth_info[1]):
                        fastq_read.append(read_nuc.lower())
                        location_meth_info[3] += 1
                    else:
                        fastq_read.append(read_nuc)
                        location_meth_info[4] += 1
                else:
                    fastq_read.append(read_nuc)
            elif read_nuc != '-':
                fastq_read.append(read_nuc)
            genome_position += 1
        if aln_profile['reference_strand'] == '-':
            fastq_read = ''.join(fastq_read[::-1])
        else:
            fastq_read = ''.join(fastq_read)
        fastq_line[1] = fastq_read
        return fastq_line

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

    def parse_aln_line(self, aln_line: List) -> Dict:
        """ Take list of line composing aln file line and returns formatted objects. Also adjust genome position
        based on strand information.
        Keyword Arguments:
            aln_line (list): list of strings in aln file format
        Returns:
            (dict) formatted aln line"""
        aln_head = aln_line[0].split('\t')
        read_index = int(aln_head[2])
        read_contig = aln_head[0].replace('>', '')
        reference_strand = aln_head[3]
        reference_sequence = aln_line[1]
        read_sequence = aln_line[2]
        if reference_strand == '-':
            reference_sequence_len = len(reference_sequence.replace('-', ''))
            read_sequence = read_sequence[::-1]
            reference_sequence = reference_sequence[::-1]
            read_index = self.reference_contig_size[read_contig] - read_index - reference_sequence_len
        return dict(read_index=read_index, read_contig=read_contig, reference_strand=reference_strand,
                    reference_sequence=reference_sequence, read_sequence=read_sequence)

    @staticmethod
    def write_fastq(output_object, fastq_line):
        for line in fastq_line:
            output_object.write(f'{line}\n')

    @property
    def get_output_objects(self):
        output_list = [open(f'{self.output_path}_meth_1.fastq', 'w')]
        if self.paired_end:
            output_list.append(open(f'{self.output_path}_meth_2.fastq', 'w'))
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
