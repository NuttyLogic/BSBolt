import random
import numpy as np
from BSB.BSB_Utils.FastaIterator import OpenFasta
from BSB.BSB_Utils.CGmapIterator import OpenCGmap
from BSB.BSB_Simulate.SimulationOutput import SimulationOutput


class SetCytosineMethylation:
    """Set cytosine methylation at CpG and CH sites. Methylation value are assigned for at + and - cytosines. After
    setting a methylation value, a random roll is performed for each read. If the roll is less than the methylation
    value the cytosine is set as methylation and if it is greater than the cytosine is unmethylated. This introduces
    some random noise into the simulation process. The number of methylated and unmethylated reads is tracked in the
    final simulation output.
    Keyword Arguments:
        reference_file (str): path to reference genome
        methylation_reference_output (str): path to output folder for simulation data
        methylation_reference (str): CGmap file with methylation profile
        methylation_profile (str): path to previously generated methylation reference
    Attributes:
        self.reference_file (str): path to reference genome
        self.methylation_reference_output (str): path to output folder for simulation data
        self.methylation_reference (str): CGmap file with methylation profile
        self.methylation_profile (str): path to previously generated methylation reference
        self.reference_dict (dict): contig sequences
        self.reference_contig_size (dict): lens of all contigs
        self.methylation_profile (dict): dict with cytosine methylation values
        self.simulation_out (dict): objects for loading and outputting simulation data
        self.cpg_distribution (np.array): beta binomial distribution for CpG sites
        self.ch_distribution (np.array): beat binomial distribution for CH sites
            """

    def __init__(self, reference_file=None, methylation_reference_output=None,
                 methylation_reference=None, methylation_profile=None):
        self.reference_file = reference_file
        self.methylation_reference_output = methylation_reference_output
        self.methylation_reference = methylation_reference
        self.reference_dict = self.get_reference_contigs
        self.reference_contig_size = {contig: len(sequence) for contig, sequence in self.reference_dict.items()}
        self.methylation_profile = methylation_profile
        self.simulation_out = None
        self.cpg_distribution = np.random.beta(.5, .5, size=5000)
        self.ch_distribution = np.random.beta(.01, .05, size=5000)

    def set_simulated_methylation(self):
        if self.methylation_reference:
            self.simulation_out = SimulationOutput(simulation_output=self.methylation_reference_output,
                                                   reference_contigs=list(self.reference_contig_size.keys()),
                                                   bsb_methylation_reference=self.methylation_reference)
            self.simulation_out.generate_simulation_directory()
        elif self.methylation_profile:
            self.simulation_out = SimulationOutput(simulation_output=self.methylation_reference_output,
                                                   reference_contigs=list(self.reference_contig_size.keys()))
            self.simulation_out.generate_simulation_directory()
            self.simulation_out.output_contig_key()
            self.get_methylation_profile_from_cgmap()
        else:
            self.simulation_out = SimulationOutput(simulation_output=self.methylation_reference_output,
                                                   reference_contigs=list(self.reference_contig_size.keys()))
            self.simulation_out.generate_simulation_directory()
            self.simulation_out.output_contig_key()
            self.set_random_cytosine_methylation()
        return self.reference_contig_size, self.simulation_out

    def get_methylation_profile_from_cgmap(self):
        current_chromosome = None
        cytosine_dict = {'Watson': {}, 'Crick': {}}
        for line in OpenCGmap(self.methylation_profile):
            chrom, nucleotide, pos, context, methlevel = line[0], line[1], int(line[2]) - 1, line[4], float(line[7])
            if not current_chromosome:
                current_chromosome = chrom
            if current_chromosome != chrom:
                self.simulation_out.output_contig_methylation_reference(cytosine_dict, chrom)
                cytosine_dict = {'Watson': {}, 'Crick': {}}
            meth_profile = dict(nucleotide=nucleotide, methylation_level=methlevel, context=context,
                                methylated_reads=0, unmethylated_reads=0)
            if nucleotide == 'G':
                cytosine_dict['Crick'][f'{chrom}:{pos}'] = meth_profile
            elif nucleotide == 'C':
                cytosine_dict['watson'][f'{chrom}:{pos}'] = meth_profile
        self.simulation_out.output_contig_methylation_reference(cytosine_dict, current_chromosome)

    @property
    def pick_cpg_methylation(self):
        """Sample from CpG distribution"""
        return random.choice(self.cpg_distribution)

    @property
    def pick_ch_methylation(self):
        """Sample from CH distribution"""
        return random.choice(self.ch_distribution)

    @staticmethod
    def random_roll(proportion_positive: float = 0.1) -> bool:
        """Return true is random < proportion_positive, used to set individual read methylation"""
        return random.random() < proportion_positive

    @property
    def get_reference_contigs(self):
        """Retrieve contig sequence from reference fasta rather than index to ensure simulation will work without
        a generated index. """
        # store contig_id and matching sequencing
        contig_dict = {}
        contig_id = None
        contig_sequence = []
        for contig_label, line in OpenFasta(self.reference_file):
            if contig_label:
                # join sequence and add to dict if new contig encounted
                if contig_id:
                    contig_dict[contig_id] = ''.join(contig_sequence)
                    contig_sequence = []
                contig_id = line.replace('>', '').split()[0]
            else:
                contig_sequence.append(line.upper())
        # join contig and add to dict after iteration finishes
        contig_dict[contig_id] = ''.join(contig_sequence)
        return contig_dict

    def set_random_cytosine_methylation(self):
        """Iterate through reference sequence, upon encountering a Cyotsine (Watson) or Guanine (Crick):
            1. get nucleotide context
            2. pull random methylation value from appropriate distribution and store value
        """
        # iterate through reference sequence
        for contig, sequence in self.reference_dict.items():
            cytosine_dict = {'Watson': {}, 'Crick': {}}
            for index, nucleotide in enumerate(sequence):
                # retrieve nucleotide context
                context = sequence[index - 1: index + 2]
                # context won't exist at ends of reference sequence
                if context:
                    # set watson methylation
                    if nucleotide == 'C':
                        c_context = context[1:]
                        methylation_profile: dict = self.get_methylation_level(c_context, nucleotide)
                        cytosine_dict['Watson'][f'{contig}:{index}'] = methylation_profile
                    # set crick methylation
                    elif nucleotide == 'G':
                        g_context = context[0:2]
                        methylation_profile: dict = self.get_methylation_level(g_context, nucleotide)
                        cytosine_dict['Crick'][f'{contig}:{index}'] = methylation_profile
            self.simulation_out.output_contig_methylation_reference(cytosine_dict, contig)

    def get_methylation_level(self, context, nucleotide):
        """Retrieve methylation level from distribution based on nucleotide context
        Arguments:
            context (str): 2BP nucleotide context
            nucleotide (str): reference nucleotide
        Return:
            methylation_profile (dict): context, methylation level, unpopulated methylation read information
            """
        if nucleotide == 'G' and context == 'GC':
            methylation_level = self.pick_cpg_methylation
        elif nucleotide == 'G' and context == 'CG':
            methylation_level = self.pick_cpg_methylation
        else:
            methylation_level = self.pick_ch_methylation
        return dict(nucleotide=nucleotide, methylation_level=methylation_level, context=context,
                    methylated_reads=0, unmethylated_reads=0)
