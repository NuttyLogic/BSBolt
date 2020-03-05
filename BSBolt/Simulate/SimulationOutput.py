import pickle
import os
import numpy as np


class SimulationOutput:

    def __init__(self, simulation_output=None, reference_contigs=None, bsb_methylation_reference=None):
        self.simulation_output = simulation_output
        self.reference_contigs = reference_contigs
        self.bsb_methylation_reference = bsb_methylation_reference
        if not self.bsb_methylation_reference:
            self.bsb_methylation_reference = self.simulation_output

    def generate_simulation_directory(self):
        simulation_directory = '/'.join(self.simulation_output.split('/')[0:-1])
        if simulation_directory:
            if not os.path.isdir(simulation_directory):
                os.makedirs(simulation_directory, exist_ok=False)

    def output_contig_methylation_reference(self, contig_key, contig_values, contig_id):
        with open(f'{self.simulation_output}.{contig_id}.pkl', 'wb') as contig_out:
            pickle.dump(contig_key, contig_out)
        np.save(f'{self.simulation_output}.{contig_id}.values.npy', contig_values)

    def load_contig(self, contig_id):
        try:
            with open(f'{self.bsb_methylation_reference}.{contig_id}.pkl', 'rb') as contig_out:
                contig_key = pickle.load(contig_out)
            contig_values = np.load(f'{self.simulation_output}.{contig_id}.values.npy')
        except FileNotFoundError:
            print(f'{contig_id} Methylation Reference not found\n')
            return None, None
        else:
            return contig_key, contig_values
