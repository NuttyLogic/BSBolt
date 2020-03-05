import random
from typing import Dict
import numpy as np
from BSBolt.Utils.FastaIterator import OpenFasta
from BSBolt.Utils.CGmapIterator import OpenCGmap
from BSBolt.Simulate.SimulationOutput import SimulationOutput


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

    def __init__(self, reference_file: str = None, methylation_reference_output: str = None,
                 methylation_reference: str = None, methylation_profile: str = None,
                 collect_ch_sites: bool = True):
        self.reference_file = reference_file
        self.methylation_reference_output = methylation_reference_output
        self.methylation_reference = methylation_reference
        self.reference_dict = self.get_reference_contigs
        self.methylation_profile = methylation_profile
        self.simulation_out = None
        self.collect_ch_sites = collect_ch_sites
        self.cpg_distribution = np.random.beta(.5, .5, size=5000)
        self.ch_distribution = np.random.beta(.01, .05, size=5000)

    def set_simulated_methylation(self):
        if self.methylation_reference:
            self.simulation_out = SimulationOutput(simulation_output=self.methylation_reference_output,
                                                   reference_contigs=list(self.reference_dict.keys()),
                                                   bsb_methylation_reference=self.methylation_reference)
            self.simulation_out.generate_simulation_directory()
        elif self.methylation_profile:
            self.simulation_out = SimulationOutput(simulation_output=self.methylation_reference_output,
                                                   reference_contigs=list(self.reference_dict.keys()))
            self.simulation_out.generate_simulation_directory()
            self.get_methylation_profile_from_cgmap()
        else:
            self.simulation_out = SimulationOutput(simulation_output=self.methylation_reference_output,
                                                   reference_contigs=list(self.reference_dict.keys()))
            self.simulation_out.generate_simulation_directory()
            self.set_random_cytosine_methylation()
        return self.reference_dict, self.simulation_out

    def get_methylation_profile_from_cgmap(self):
        current_chromosome = None
        meth_keys = {}
        cytosine_values = {}
        cytosine_positions = {}
        for line in OpenCGmap(self.methylation_profile):
            chrom, nucleotide, pos, context, methlevel = line[0], line[1], int(line[2]) - 1, line[4], float(line[7])
            context = 1 if context == 'CG' else 0
            if not self.collect_ch_sites and not context:
                continue
            nucleotide = 1 if nucleotide == 'G' else 0
            # nucleotide, methylation_level, context, methylated_reads, unmethylated_reads
            meth_profile = np.assarray([nucleotide, methlevel, context, 0, 0, 1])
            profile_id = f'{chrom}:{pos}'
            if chrom not in meth_keys:
                meth_keys[chrom] = {profile_id: 0}
                cytosine_values[chrom] = [meth_profile]
                cytosine_positions[chrom] = 1
            else:
                meth_keys[chrom][profile_id] = cytosine_positions[chrom]
                cytosine_values[chrom].append(meth_profile)
                cytosine_positions[chrom] += 1
        for contig, meth_key in meth_keys.items():
            self.simulation_out.output_contig_methylation_reference(meth_key, np.array(cytosine_values[contig]), contig)

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
    def get_reference_contigs(self) -> Dict[str, str]:
        """Retrieve contig sequence from reference fasta rather than index to ensure simulation will work without
        a generated index. """
        # store contig_id and matching sequencing
        contig_dict = {}
        contig_id = None
        contig_sequence = []
        for contig_label, line in OpenFasta(self.reference_file):
            if contig_label:
                # join sequence and add to dict if new contig encountered
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
        """Iterate through reference sequence, upon encountering a Cytosine (Watson) or Guanine (Crick):
            1. get nucleotide context
            2. pull random methylation value from appropriate distribution and store value
        """
        # iterate through reference sequence
        for contig, sequence in self.reference_dict.items():
            cytosine_key = {}
            position_count = 0
            cytosine_values = []
            for index, nucleotide in enumerate(sequence):
                # retrieve nucleotide context
                context = sequence[index - 1: index + 2]
                # context won't exist at ends of reference sequence
                if context:
                    # set watson methylation
                    methylation_profile = []
                    if nucleotide == 'C':
                        context = 1 if context[1:] == 'CG' else 0
                        methylation_profile: list = self.get_methylation_level(context, 0)
                    # set crick methylation
                    elif nucleotide == 'G':
                        context = 1 if context[0:2] == 'CG' else 0
                        methylation_profile: list = self.get_methylation_level(context, 1)
                    if methylation_profile:
                        # if no CH information don't add info
                        if not self.collect_ch_sites and not methylation_profile[2]:
                            continue
                        cytosine_key[f'{contig}:{index}'] = position_count
                        cytosine_values.append(np.array(methylation_profile))
                        position_count += 1
            self.simulation_out.output_contig_methylation_reference(cytosine_key, np.array(cytosine_values), contig)

    def get_methylation_level(self, context, nucleotide) -> list:
        """Retrieve methylation level from distribution based on nucleotide context
        Arguments:
            context (int): 2BP nucleotide context
            nucleotide (int): reference nucleotide
        Return:
            methylation_profile (np.ndarray): context, methylation level, unpopulated methylation read information
            """
        if context:
            methylation_level = self.pick_cpg_methylation
        else:
            methylation_level = self.pick_ch_methylation
        return [nucleotide, methylation_level, context, 0, 0, 0]
