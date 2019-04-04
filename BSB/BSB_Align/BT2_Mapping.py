import multiprocessing
from BSB.BSB_Align.AlignmentHelpers import launch_bowtie2_mapping


class Bowtie2Alignment:
    """Class to launch bowtie2 alingment,
        Keyword Arguments:
            fastq1 (str): path to fastq file
            fastq2 (str): path to fastq1 mate pair
            undirectional_library (bool): perform undirectional alignment
            bowtie2_commands (list of str): commands to use for bowtie2 alignment
            bsb_database(str): path to bsb 2 directory
            bowtie2_path(str): path to bowtie2 executable
            output_path(str): output prefix of output file
            non_converted_output (bool): map non-converted reads
        Attributes:
            self.bowtie2_mapping (dict): argument for launching bowtie2 mapping
            self.fastq1 (str): path to fastq file
            self.fastq2 (str): path to fastq1 mate pair
            self.undirectional_library (bool): perform undirectional alignment
            self.bsb_database(str): path to bsb 2 directory
            self.output_path (str): output prefix of output file
            self.non_converted_output (bool): map non-converted reads
            self.paired_end (bool): paired end fastq files
            self.mapping_commands (list): list of mapping commands to launch
            """

    def __init__(self, fastq1=None, fastq2=None, undirectional_library=False, bowtie2_commands=None,
                 bsb_database=None, bowtie2_path=None, output_path=None, non_converted_output=False):
        self.mapping_dict = dict(bowtie2_path=bowtie2_path, fastq1=fastq1, fastq2=fastq2,
                                 bowtie2_commands=bowtie2_commands)
        self.tab_stream_commands = dict(fastq1=fastq1, fastq2=fastq2)
        assert isinstance(undirectional_library, bool)
        self.undirectional_library = undirectional_library
        self.paired_end = False
        if fastq2:
            self.paired_end = True
        self.output_path = output_path
        self.non_converted_output = non_converted_output
        self.bowtie2_commands = bowtie2_commands
        self.bsb_database = bsb_database
        self.mapping_commands = self.get_mapping_commands

    @property
    def get_mapping_commands(self):
        """ Formats and returns mapping commands, setting correct base substitutions for mapping
        Returns:
             list of formatted mapping commands
        """
        # set default variable bases
        r_base1, r_base2 = 'C', 'G'
        if self.non_converted_output:
            r_base1, r_base2 = None, None
        # Launch W_C2T mapping command
        W_C2T_mapping = dict(self.mapping_dict)
        W_C2T_mapping.update(dict(replacement_base1=r_base1, replacement_base2=r_base2,
                                  bowtie2_database=f'{self.bsb_database}W_C2T'))
        # Launch C_C2T mapping command
        C_C2T_mapping = dict(self.mapping_dict)
        C_C2T_mapping.update(dict(replacement_base1=r_base1, replacement_base2=r_base2,
                                  bowtie2_database=f'{self.bsb_database}C_C2T'))
        if self.undirectional_library:
            # Launch W_G2A mapping command
            W_G2A_mapping = dict(self.mapping_dict)
            W_G2A_mapping.update(dict(replacement_base1=r_base2, replacement_base2=r_base1,
                                      bowtie2_database=f'{self.bsb_database}W_G2A'))
            # Launch C_G2A mapping command
            C_G2A_mapping = dict(self.mapping_dict)
            C_G2A_mapping.update(dict(replacement_base1=r_base2, replacement_base2=r_base1,
                                      bowtie2_database=f'{self.bsb_database}C_G2A'))
            return [(W_C2T_mapping, 'W_C2T'), (C_C2T_mapping, 'C_C2T'),
                    (W_G2A_mapping, 'W_G2A'), (C_G2A_mapping, 'C_G2A')]
        return [(W_C2T_mapping, 'W_C2T'), (C_C2T_mapping, 'C_C2T')]

    def launch_bt2_mapping(self):
        """Launch asynchronous mapping commands, mapping commands will be executed in parallel and block until
        all processes are complete"""
        # pool thread limit of four, can't go above this, can limit if resources are limited
        pool_threads = 4
        # start pool
        pool = multiprocessing.Pool(processes=pool_threads)
        for mapping_command, genome_database_label in self.mapping_commands:
            pool.apply_async(launch_bowtie2_mapping,
                             kwds={'bowtie2_stream_kwargs': mapping_command,
                                   'output_path': self.output_path,
                                   'genome_database_label': genome_database_label},
                             error_callback=self.propagate_error)
        pool.close()
        # join and wait for processes to complete
        pool.join()

    def propagate_error(self, error):
        raise error
