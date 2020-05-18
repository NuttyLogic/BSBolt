import random
from typing import Dict, List
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

    def __init__(self, reference_file: str = None, sim_dir: str = None, methylation_reference: str = None,
                 cgmap: str = None, collect_ch_sites: bool = True, overwrite_db: bool = False):
        self.reference = self.get_reference(reference_file)
        self.reference_file = reference_file
        self.sim_dir = sim_dir
        self.sim_db: SimulationOutput = None
        self.initialize_methylation_reference(methylation_reference, cgmap)
        self.methylation_reference = methylation_reference
        self.collect_ch_sites = collect_ch_sites
        self.overwrite_db = overwrite_db
        self.cpg_distribution = np.random.beta(.5, .5, size=5000)
        self.ch_distribution = np.random.beta(.01, .05, size=5000)

    def initialize_methylation_reference(self, methylation_reference, cgmap):
        """If profile exists then assume already run"""
        if methylation_reference:
            meth_in = SimulationOutput(sim_dir=methylation_reference)
            meth_out = SimulationOutput(sim_dir=self.sim_dir)
            meth_out.generate_simulation_directory()
            for contig in self.reference.keys():
                for contig_key, contig_values in meth_in.load_contig(contig):
                    if contig_key:
                        meth_out.output_contig(contig_key, contig_values, contig)
            self.sim_db = meth_out
        elif cgmap:
            self.sim_db = SimulationOutput(sim_dir=self.sim_dir)
            self.sim_db.generate_simulation_directory()
            self.sim_db_from_cgmap(cgmap)
        else:
            self.sim_db = SimulationOutput(sim_dir=self.sim_dir)
            self.sim_db.generate_simulation_directory()

    def sim_db_from_cgmap(self, cgmap: str):
        meth_keys = {}
        cytosine_values = {}
        cytosine_positions = {}
        for line in OpenCGmap(cgmap):
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
            self.sim_db.output_contig(meth_key, np.array(cytosine_values[contig]), contig)

    @property
    def pick_cpg_methylation(self):
        """Sample from CpG distribution"""
        return random.choice(self.cpg_distribution)

    @property
    def pick_ch_methylation(self):
        """Sample from CH distribution"""
        return random.choice(self.ch_distribution)

    @staticmethod
    def get_reference(reference_file: str) -> Dict[str, str]:
        """Retrieve contig sequence from reference fasta rather than index to ensure simulation will work without
        a generated index. """
        # store contig_id and matching sequencing
        reference = {}
        contig_id = None
        contig_sequence = []
        for contig_label, line in OpenFasta(reference_file):
            if contig_label:
                # join sequence and add to dict if new contig encountered
                if contig_id:
                    reference[contig_id] = ''.join(contig_sequence)
                    contig_sequence = []
                contig_id = line.replace('>', '').split()[0]
            else:
                contig_sequence.append(line.upper())
        # join contig and add to dict after iteration finishes
        reference[contig_id] = ''.join(contig_sequence)
        return reference

    def get_contig_methylation(self, contig):
        if self.overwrite_db:
            contig_profile = None
        else:
            contig_profile = self.sim_db.load_contig(contig)
        if contig_profile:
            return contig_profile
        else:
            return self.set_random_cytosine_methylation(contig)

    def set_variant_methylation(self, sim_data, contig_profile, current_contig):
        """Variants are always random, so set random methylation"""
        ref_seq = self.reference[current_contig]
        for pos, variant_info in sim_data.items():
            # don't set methylation for last base
            if pos + 2 > len(ref_seq):
                continue
            if variant_info['indel'] == -1:
                continue
            elif variant_info['indel'] == 1:
                insert_context = f'{ref_seq[pos]}{variant_info["alt"]}{ref_seq[pos + 1]}'
                for insert_pos, nucleotide in enumerate(variant_info['alt']):
                    local_context = insert_context[insert_pos: insert_pos + 3]
                    meth_profile = self.get_random_methylation_profile(nucleotide, local_context)
                    if meth_profile:
                        contig_profile[f'{pos}_+_{insert_pos + 1}'] = meth_profile
            else:
                for nucleotide in variant_info['iupac']:
                    if nucleotide == variant_info['reference']:
                        continue
                    else:
                        context = f'{ref_seq[pos - 1]}{nucleotide}{ref_seq[pos + 1]}'
                        meth_profile = self.get_random_methylation_profile(nucleotide, context)
                        if meth_profile:
                            contig_profile[f'{pos}_{variant_info["reference"]}_{nucleotide}'] = meth_profile

    def set_random_cytosine_methylation(self, contig):
        """Iterate through reference sequence, upon encountering a Cytosine (Watson) or Guanine (Crick):
            1. get nucleotide context
            2. pull random methylation value from appropriate distribution and store value
        """
        # iterate through reference sequence
        sequence = self.reference[contig]
        cytosine_profile = {}
        for pos, nucleotide in enumerate(sequence):
            # retrieve nucleotide context
            context = sequence[pos - 1: pos + 2]
            # context won't exist at ends of reference sequence
            if context:
                # set watson methylation
                methylation_profile = self.get_random_methylation_profile(nucleotide, context)
                if methylation_profile:
                    # if no CH information don't add info
                    if not self.collect_ch_sites and not methylation_profile[2]:
                        continue
                    cytosine_profile[f'{pos}'] = tuple(methylation_profile)
        return cytosine_profile

    def get_random_methylation_profile(self, nucleotide, context):
        if nucleotide == 'C':
            n_context = 1 if context[1:] == 'CG' else 0
            return self.get_methylation_level(n_context, 0)
        # set crick methylation
        elif nucleotide == 'G':
            n_context = 1 if context[0:2] == 'CG' else 0
            return self.get_methylation_level(n_context, 1)
        return None

    def get_methylation_level(self, context, nucleotide) -> tuple:
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
        return nucleotide, methylation_level, context
