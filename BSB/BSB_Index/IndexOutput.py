import os
import subprocess
import pickle
import gzip

from BSB.BSB_Utils.UtilityFunctions import reverse_complement


class IndexOutput:
    """Class to output processed sequence index sequence and launch external bowtie2 build commands.
    Keyword Arguments
        genome_database (str):
        bowtie2_path (str):
        bowtie2_threads (int):
    Attributes:
        self.genome_database (str): formatted str to ensure proper file output
        self.bowtie2_path (str): bowtie2 command if in path or path to executable
        self.bowtie2_threads (int): thread count for bowtie2
        self.database_output_object (Dict[str, TextIO]): Contains TextIO output objects for writing processed DNA
                                                         sequence
    """

    def __init__(self, genome_database=None, bowtie2_path=None, bowtie2_threads=1):
        assert isinstance(genome_database, str), 'Genome Database Path Invalid, Must be a String'
        assert isinstance(bowtie2_path, str), 'Bowtie2 Path Invalid, Must be a String'
        assert isinstance(bowtie2_threads, int), 'Bowtie2 Threads Invalid, Must be Integer'
        # format genome_database path
        self.genome_database = self.generate_genome_directory(genome_database)
        self.bowtie2 = bowtie2_path
        self.bowtie2_threads = str(bowtie2_threads)
        # set output object
        self.database_output = open(f'{self.genome_database}BSB_ref.fa', 'w')

    @staticmethod
    def generate_genome_directory(genome_database):
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

    def write_contig_sequence(self, contig_id, contig_sequence):
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

    def build_bowtie2_index(self):
        """Launch external bowtie2 commands for 4 processed reference files and collect external stdout for log file
        """
        # format output and input
        index_input = f'{self.genome_database}BSB_ref.fa'
        index_output = f'{self.genome_database}BSB_ref'
        # collect external command
        bowtie_command = [f'{self.bowtie2}-build',
                          '--threads', self.bowtie2_threads,
                          '-f', index_input,
                          index_output]
        # open file to collect stdout
        bowtie2_index_log = open(f'{self.genome_database}.bt2_index.log', 'w')
        # run external command
        subprocess.run(args=bowtie_command, stdout=bowtie2_index_log)

    def output_contig_sequence(self, contig_id, contig_sequence):
        """Outputs serialized version of contig sequence
            Arguments:
                contig_id (str): contig label
                contig_sequence (object): str of DNA sequence
        """
        with open(f'{self.genome_database}{contig_id}.pkl', 'wb') as contig:
            return pickle.dump(contig_sequence, contig)

    def output_mappable_regions(self, mappable_regions):
        """Outputs mappable regions
        Arguments:
            mappable_regions (list): list of bed formatted strings
                """
        with gzip.open(f'{self.genome_database}mappable_regions.txt.gz', 'wb') as mappable_regions_output:
            for line in mappable_regions:
                mappable_regions_output.write(line.encode('UTF-8'))
