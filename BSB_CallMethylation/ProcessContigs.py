import multiprocessing
import gzip
import pysam
import time
from tqdm import tqdm
from BSB_CallMethylation.CallMethylation import CallMethylation


class MethylationCallingError(Exception):
    """Error in methylation calling process"""
    pass


def call_contig_methylation(return_dict, call_methylation_kwargs):
    """ Wrapper to initialize methylation calling
    Arguments:
        return_dict (multiprocessing.mangager.dict): Dictionary to return processed methylation calls
        call_methylation_kwargs (dict): dict of argument for CallMethylation class
    """
    contig_methylation_call = CallMethylation(**call_methylation_kwargs)
    contig_methylation_call.call_methylation()
    assert isinstance(contig_methylation_call, CallMethylation)
    return_dict[call_methylation_kwargs['contig']] = 1


class ProcessContigs:
    """
    Keyword Arguments:
        input_file (str): str to input bam/sam file
        genome_database (str): str to genome directory
        output_prefix (str): output prefix for ATCGmap, CGmap, and WIG files
        remove_sx_reads (bool): remove incompletely converted reads
        ignore_overlap (bool):  ignore overlapping reads
        text_output (bool): output compressed or plain text files
        remove_ccgg (bool): don't call CCGG sequences
        min_read_depth (int): default = 1
        max_read_depth (int): default = 8000
        threads (int): 1, if one watcher and processing on same thread else separated
        min_base_quality (int): minimum base quality for base to be considered
    Attributes:
        self.input_file (str): path to input bam/sam file
        self.input_bam (pysam.Samfile): pysam object to retrieve pileup information from .bam file
        self.text_output (bool): plain text output
        self.output_prefix (str): output prefix
        self.threads (int): , if one thread watcher and processing on same thread else separated
        self.call_methylation_kwargs (dict): dict of kwargs to pass to CallMethylation class
        self.min_read_depth (int): default = 1, minimum read depth to output values
        self.contigs (list): list of contigs in self.input_bam
        self.return_dict (multiprocessing.manager.dict): thread safer dictionary
        self.output_objects (Dict[str, TextIO]): dict of output objects
    """

    def __init__(self, input_file=None, genome_database=None, output_prefix=None,
                 remove_sx_reads=True, ignore_overlap=False, text_output=False, remove_ccgg=False,
                 min_read_depth=10, max_read_depth=8000, threads=1, verbose=True, min_base_quality=0):
        assert isinstance(input_file, str), 'Path to input file not valid'
        assert isinstance(text_output, bool), 'Not valid bool'
        assert isinstance(threads, int), 'Threads must be specified with integer'
        if output_prefix:
            assert isinstance(output_prefix, str)
        self.input_file = input_file
        self.input_bam = pysam.Samfile(input_file, 'rb')
        self.text_output = text_output
        self.output_prefix = output_prefix
        self.threads = threads
        self.call_methylation_kwargs = dict(input_file=input_file,
                                            genome_database=genome_database,
                                            remove_sx_reads=remove_sx_reads,
                                            ignore_overlap=ignore_overlap,
                                            remove_ccgg=remove_ccgg,
                                            min_read_depth=min_read_depth,
                                            max_read_depth=max_read_depth,
                                            min_base_quality=min_base_quality)
        self.min_read_depth = min_read_depth
        self.methylation_calling = True
        self.contigs = self.get_contigs
        self.return_dict = None
        self.return_list = None
        self.pool = None
        self.verbose = verbose
        self.output_objects = self.get_output_objects
        self.process_contigs()
        self.watch_pool()

    @property
    def get_contigs(self):
        """
        get list of contigs in input file, threads set across contigs
        """
        return [contig for contig in self.input_bam.references]

    def process_contigs(self):
        """Launches a processing pool to call methylation values across the input file contigs
        """
        # initialize manager
        manager = multiprocessing.Manager()
        # get return dictionary
        self.return_list = manager.list()
        self.return_dict = manager.dict()
        # threads for methylation calling, if one thread use thread for calling and watching
        pool_threads = self.threads - 1 if self.threads != 1 else 1
        # start pool
        self.pool = multiprocessing.Pool(processes=pool_threads)
        # for contig call methylation and return values to dict
        for contig in self.contigs:
            contig_kwargs = dict(self.call_methylation_kwargs)
            contig_kwargs.update(dict(contig=contig, return_list=self.return_list))
            self.pool.apply_async(call_contig_methylation,
                                  args=[self.return_dict, contig_kwargs],
                                  error_callback=self.methylation_process_error)
        self.pool.close()

    def methylation_process_error(self, error):
        """Raise if exception thrown in methylation calling process"""
        self.methylation_calling = False
        raise MethylationCallingError(error)

    def watch_pool(self):
        """Watch self.return_dict and process return methylation values. Contigs are processed in order so buffer can
        become large if first contig is large, ie. Human Chr1
        """
        completed_contigs = set()
        if self.verbose:
            pbar = tqdm(total=len(self.contigs), desc='Processing Contigs')
        while self.methylation_calling:
            try:
                if self.return_list:
                    methylation_lines: list = self.return_list.pop(0)
            except IndexError:
                # if contig is missing sleep
                if len(self.return_dict) == len(self.contigs) and not self.return_list:
                    self.methylation_calling = False
                else:
                    for contig in self.return_dict:
                        if contig not in completed_contigs:
                            if self.verbose:
                                pbar.update(1)
                            completed_contigs.add(contig)
                time.sleep(2)
            else:
                # write output
                self.write_output(methylation_lines, self.contigs[0])
        if self.verbose:
            pbar.close()

    def write_output(self, methylation_lines, contig):
        """Give a list of methylation call dicts, output formatted line
        Arguments:
            methylation_lines (list): list of dict containing methylation call information
            contig (str): contig label
        """
        # write wig contig designation
        self.write_line(self.output_objects['wig'], f'variableStep chrom={contig}\n')
        for meth_line in methylation_lines:
            # write ATCGmap line
            self.write_line(self.output_objects['ATCGmap'], self.format_atcg(meth_line))
            # if methylation level greater than or equal to min_read_depth output CGmap and wig lines
            if meth_line['all_cytosines'] >= self.min_read_depth:
                self.write_line(self.output_objects['CGmap'], self.format_cgmap(meth_line))
                self.write_line(self.output_objects['wig'], self.format_wig(meth_line))

    @staticmethod
    def format_atcg(meth_line):
        """ATCGmap line formatting
        """
        ATCGmap_line = f'{meth_line["chrom"]}\t{meth_line["nucleotide"]}\t{meth_line["pos"]}\t{meth_line["context"]}' \
                       f'\t{meth_line["subcontext"]}\t{meth_line["forward_counts"]}\t{meth_line["reverse_counts"]}' \
                       f'\t{meth_line["meth_level"]}\n'
        return ATCGmap_line

    @staticmethod
    def format_cgmap(meth_line):
        """CGmap line formatting
        """
        CGmap_line = f'{meth_line["chrom"]}\t{meth_line["nucleotide"]}\t{meth_line["pos"]}' \
                     f'\t{meth_line["context"]}\t{meth_line["subcontext"]}\t{meth_line["meth_level"]}' \
                     f'\t{meth_line["meth_cytosines"]}\t{meth_line["all_cytosines"]}\n'
        return CGmap_line

    @staticmethod
    def format_wig(meth_line):
        """WIG line formatting
        """
        wiggle_line = f'{meth_line["pos"]}\t{meth_line["meth_level"]:.2f}\n'
        if meth_line["nucleotide"] == 'C':
            wiggle_line = f'{meth_line["pos"]}\t{-meth_line["meth_level"]:.2f}\n'
        return wiggle_line

    @property
    def get_output_objects(self):
        output_objects = {}
        output_path = self.input_file
        if self.output_prefix:
            output_path = self.output_prefix
        if self.text_output:
            output_objects['ATCGmap'] = open(f'{output_path}.ATCGmap', 'w')
            output_objects['CGmap'] = open(f'{output_path}.CGmap', 'w')
            output_objects['wig'] = open(f'{output_path}.wig', 'w')
            output_objects['wig'].write('track type=wiggle_o\n')
        else:
            output_objects['ATCGmap'] = gzip.open(f'{output_path}.ATCGmap.gz', 'wb')
            output_objects['CGmap'] = gzip.open(f'{output_path}.CGmap.gz', 'wb')
            output_objects['wig'] = gzip.open(f'{output_path}.wig.gz', 'wb')
            output_objects['wig'].write('track type=wiggle_o\n'.encode('utf-8'))
        return output_objects

    def write_line(self, output_object, line):
        """ Outputs line, and encodes before output if necessary
        Arguments
            output_object (TextIO/GZipIO): output object
            line: formatted line to write
        """
        if self.text_output:
            output_object.write(line)
        else:
            output_object.write(line.encode('utf-8'))
