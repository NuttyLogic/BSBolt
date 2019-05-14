import pickle
import os


class SimulationOutput:

    def __init__(self, simulation_output=None, reference_contigs=None, bsb_methylation_reference=None):
        self.simulation_output = simulation_output
        self.reference_contigs = reference_contigs
        self.bsb_methylation_reference = bsb_methylation_reference
        if not self.bsb_methylation_reference:
            self.bsb_methylation_reference = self.simulation_output

    def generate_simulation_directory(self):
        simulation_directory = '/'.join(self.simulation_output.split('/')[0:-1])
        if not os.path.isdir(simulation_directory):
            os.makedirs(simulation_directory, exist_ok=False)

    def output_contig_key(self):
        with open(f'{self.simulation_output}.genome_index.pkl', 'wb') as ref_out:
            pickle.dump(self.reference_contigs, ref_out)

    def output_contig_methylation_reference(self, contig_methylation, contig_id):
        with open(f'{self.simulation_output}.{contig_id}.pkl', 'wb') as contig_out:
            pickle.dump(contig_methylation, contig_out)

    def load_contig(self, contig_id):
        with open(f'{self.bsb_methylation_reference}.{contig_id}.pkl', 'rb') as contig_out:
            contig_methylation = pickle.load(contig_out)
            return contig_methylation

    def load_contig_key(self):
        with open(f'{self.bsb_methylation_reference}.genome_index.pkl', 'rb') as ref_out:
            contig_keys = pickle.load(ref_out)
            return contig_keys
