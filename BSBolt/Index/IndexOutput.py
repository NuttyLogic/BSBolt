import os
import subprocess
import pickle
import gzip
from typing import Dict, List, Union

from BSBolt.Utils.UtilityFunctions import get_external_paths, reverse_complement


class IndexOutput:
    """Class to output processed sequence index sequence and launch external index commands.

    """

    def __init__(self, genome_database: str = None):
        # format genome_database path
        bwa_path, _ = get_external_paths()
        self.genome_database = self.generate_genome_directory(genome_database)
        self.bwa_path = bwa_path
        # set output object
        self.database_output = open(f'{self.genome_database}BSB_ref.fa', 'w')

    @staticmethod
    def generate_genome_directory(genome_database: str) -> str:
        """ Make directory if it doesn't exist, add / to output to ensure proper formatting
        Arguments:
            genome_database (str): output folder
        Returns:
             genome_database (str): formatted output path
        """
        if not os.path.isdir(genome_database):
            os.makedirs(genome_database, exist_ok=False)
        if not genome_database.endswith('/'):
            genome_database = f'{genome_database}/'
        return genome_database

    def write_contig_sequence(self, contig_id: str, contig_sequence: str):
        """ Writes formatted DNA sequence. Write possible outputs for Watson and Crick strands.
        Arguments:
            contig_id (str): contig label
            contig_sequence (str): str of DNA sequence
        """
        # write forward bisulfite converted sequence
        self.database_output.write(f'>{contig_id}\n')
        self.database_output.write(f'{contig_sequence.replace("C", "T").replace("c", "t")}\n')

        # write reverse bisulfite converted sequence
        reverse_contig_sequence = reverse_complement(contig_sequence)
        self.database_output.write(f'>{contig_id}_crick_bs\n')
        self.database_output.write(f'{reverse_contig_sequence.replace("C", "T").replace("c", "t")}\n')

    def build_index(self):
        """Launch external commands for 4 processed reference files and collect external stdout for log file
        """
        # format output and input
        ref_file = f'{self.genome_database}BSB_ref.fa'
        # collect external command
        indx_command = [f'{self.bwa_path}', 'index', ref_file]
        # run external command
        subprocess.run(args=indx_command)

    def output_contig_sequence(self, contig_id: str, contig_sequence: Union[str, Dict[str, int]]):
        """Outputs serialized version of contig sequence
            Arguments:
                contig_id (str): contig label
                contig_sequence (object): str of DNA sequence
        """
        with open(f'{self.genome_database}{contig_id}.pkl', 'wb') as contig:
            return pickle.dump(contig_sequence, contig)

    def output_mappable_regions(self, mappable_regions: List[str]):
        """Outputs mappable regions
        Arguments:
            mappable_regions (list): list of bed formatted strings
                """
        with gzip.open(f'{self.genome_database}mappable_regions.txt.gz', 'wb') as mappable_regions_output:
            for line in mappable_regions:
                mappable_regions_output.write(line.encode('UTF-8'))
