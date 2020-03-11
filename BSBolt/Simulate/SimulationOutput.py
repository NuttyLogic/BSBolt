import pickle
import os
import numpy as np


class SimulationOutput:

    def __init__(self, sim_dir=None):
        self.sim_dir = sim_dir

    def generate_simulation_directory(self):
        simulation_directory = '/'.join(self.sim_dir.split('/')[0:-1])
        if simulation_directory:
            if not os.path.isdir(simulation_directory):
                os.makedirs(simulation_directory, exist_ok=False)

    def output_contig(self, contig_key, contig_values, contig_id):
        if contig_id:
            with open(f'{self.sim_dir}.{contig_id}.pkl', 'wb') as contig_out:
                pickle.dump(contig_key, contig_out)
            np.save(f'{self.sim_dir}.{contig_id}.values.npy', contig_values)

    def load_contig(self, contig_id):
        try:
            with open(f'{self.sim_dir}.{contig_id}.pkl', 'rb') as contig_out:
                contig_key = pickle.load(contig_out)
            contig_values = np.load(f'{self.sim_dir}.{contig_id}.values.npy')
        except FileNotFoundError:
            print(f'{contig_id} Methylation Reference not found\n')
            return None, None
        else:
            return contig_key, contig_values
