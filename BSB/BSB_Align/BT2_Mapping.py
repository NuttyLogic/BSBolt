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
        r_base1, r_base2 = 'C', 'G'
        if non_converted_output:
            r_base1, r_base2 = None, None
        self.mapping_dict = dict(bowtie2_path=bowtie2_path,
                                 fastq1=fastq1, fastq2=fastq2,
                                 bowtie2_commands=bowtie2_commands,
                                 replacement_base1=r_base1, replacement_base2=r_base2,
                                 bowtie2_database=bsb_database)
        self.tab_stream_commands = dict(fastq1=fastq1, fastq2=fastq2, undirectional_library=undirectional_library)
        self.output_path = output_path

    def launch_bt2_mapping(self):
        """Launch asynchronous mapping commands, mapping commands will be executed in parallel and block until
        all processes are complete"""
        # pool thread limit of four, can't go above this, can limit if resources are limited
        pool_threads = 4
        # start pool
        try:
            pool = multiprocessing.Pool(processes=pool_threads)
        except OSError as e:
            print('Shared memory exhausted, please reduce the number of Bowtie2 Threads, or '
                  'increase the memory allocation')
            raise e
        pool.apply_async(launch_bowtie2_mapping,
                         kwds={'bowtie2_stream_kwargs': self.mapping_dict,
                               'output_path': self.output_path},
                         error_callback=self.propagate_error)
        pool.close()
        # join and wait for processes to complete
        pool.join()

    def propagate_error(self, error):
        raise error
