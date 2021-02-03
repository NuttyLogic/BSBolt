import pickle
import os


class SimulationOutput:

    def __init__(self, sim_dir=None):
        self.sim_dir = sim_dir

    def generate_simulation_directory(self):
        simulation_directory = '/'.join(self.sim_dir.split('/')[0:-1])
        if simulation_directory:
            if not os.path.isdir(simulation_directory):
                os.makedirs(simulation_directory, exist_ok=False)

    def output_contig(self, contig_profile, contig_id, values=False, variant=False):
        if contig_id:
            contig_label = contig_id
            if values:
                contig_label = f'{contig_id}_values'
            elif variant:
                contig_label = f'{contig_id}_variants'
            with open(f'{self.sim_dir}.{contig_label}.pkl', 'wb') as contig_out:
                pickle.dump(contig_profile, contig_out)

    def load_contig(self, contig_id, values=False):
        contig_label = contig_id
        if values:
            contig_label = f'{contig_id}_values'
        try:
            with open(f'{self.sim_dir}.{contig_label}.pkl', 'rb') as contig_out:
                contig_profile = pickle.load(contig_out)
        except FileNotFoundError:
            return None
        else:
            return contig_profile

