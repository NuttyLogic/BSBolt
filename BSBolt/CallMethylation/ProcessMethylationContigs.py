import io
import gzip
import multiprocessing
from typing import Dict, List, Tuple, Union
import pysam
from tqdm import tqdm
from BSBolt.CallMethylation.CallValues import CallMethylationValues


class MethylationCallingError(Exception):
    """Error in methylation calling process"""
    pass


def call_contig_methylation(completed_contigs, call_methylation_kwargs):
    """ Wrapper to initialize methylation calling
    Arguments:
        completed_contigs (multiprocessing.manager.list): List of contigs with completed methylation calls
        call_methylation_kwargs (dict): dict of argument for CallMethylation class
    """
    contig_methylation_call = CallMethylationValues(**call_methylation_kwargs)
    contig_methylation_call.call_methylation()
    assert isinstance(contig_methylation_call, CallMethylationValues)
    completed_contigs.append(call_methylation_kwargs['contig'])


class ProcessContigs:
    """
    Multi-threaded contig processing wrapper. Passes thread safe queue to workers and outputs values as CGmap file.

    Params:

   * *input_file (str)*: str to input bam/sam file
   * *genome_database (str)*: str to genome directory
   * *output_prefix (str)*: output prefix for  CGmap
   * *ignore_overlap (bool)*:  ignore overlapping reads, [True]
   * *text_output (bool)*: output compressed or plain text files, [False]
   * *remove_ccgg (bool)*: don't call CCGG sequences, [False]
   * *min_read_depth (int)*: default = minimum read depth to call methylation values, [1]
   * *max_read_depth (int)*: maximum read depth for pileup, [8000]
   * *threads (int)*: , if one watcher and processing on same thread else separated, [1]
   * *min_base_quality (int)*: minimum base quality for base to be considered, [10]
   * *min_mapping_quality (int)*: minimum mapping quality for an alignment to be considered, [10]
   * *verbose (bool)*: Verbose processing, [False]
   * *cg_only (bool)*: only return CG sites to queue, [False]
   * *ignore_oprhans (bool)*: ignore orphaned reads (not properly paired), [True]

    Usage:

    ```python
    process_values = ProcessContigs(**kwargs)
    process_values.process_contigs()
    ```
    """

    def __init__(self, input_file: str = None, genome_database: str = None, output_prefix: str = None,
                 ignore_overlap: bool = True, text_output: bool = False, remove_ccgg: bool = False,
                 min_read_depth: int = 10, max_read_depth: int = 8000, threads: int = 1, verbose: bool = True,
                 min_base_quality: int = 10, min_mapping_quality: int = 10,
                 ATCGmap: bool = False, cg_only: bool = True, ignore_orphans: bool = False):
        assert isinstance(input_file, str), 'Path to input file not valid'
        assert isinstance(text_output, bool), 'Not valid bool'
        assert isinstance(threads, int), 'Threads must be specified with integer'
        if output_prefix:
            assert isinstance(output_prefix, str)
        self.input_file = input_file
        try:
            self.input_bam = pysam.Samfile(input_file, 'rb', require_index=True)
        except IOError:
            print('Generating Index File')
            pysam.index(input_file)
            self.input_bam = pysam.Samfile(input_file, 'rb', require_index=True)
        self.text_output = text_output
        self.output_prefix = output_prefix
        self.threads = threads
        self.call_methylation_kwargs = dict(input_file=input_file,
                                            genome_database=genome_database,
                                            ignore_overlap=ignore_overlap,
                                            ignore_orphans=ignore_orphans,
                                            remove_ccgg=remove_ccgg,
                                            max_read_depth=max_read_depth,
                                            min_base_quality=min_base_quality,
                                            min_mapping_quality=min_mapping_quality,
                                            cg_only=cg_only)
        self.min_read_depth = min_read_depth
        self.ATCGmap = ATCGmap
        self.methylation_calling = True
        self.contigs = self.get_contigs
        self.completed_contigs = None
        self.return_queue = None
        self.pool = None
        self.verbose = verbose
        self.output_objects = self.get_output_objects
        self.methylation_stats = {'CG_meth': 0, 'CG_all': 0, 'CH_meth': 0, 'CH_all': 0}

    @property
    def get_contigs(self) -> List[str]:
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
        self.return_queue = manager.Queue(maxsize=20)
        self.completed_contigs = manager.list()
        # threads for methylation calling, if one thread use thread for calling and watching
        pool_threads = self.threads - 1 if self.threads != 1 else 1
        # start pool
        self.pool = multiprocessing.Pool(processes=pool_threads)
        # for contig call methylation and return values to dict
        for contig in self.contigs:
            contig_kwargs = dict(self.call_methylation_kwargs)
            contig_kwargs.update(dict(contig=contig, return_queue=self.return_queue))
            self.pool.apply_async(call_contig_methylation,
                                  args=[self.completed_contigs, contig_kwargs],
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
        contigs_complete = 0
        pbar = None
        if self.verbose:
            pbar = tqdm(total=len(self.contigs), desc='Processing Contigs')
        while self.methylation_calling:
            methylation_lines: list = self.return_queue.get(block=True)
            if self.verbose:
                if len(self.completed_contigs) != contigs_complete:
                    update_number = len(self.completed_contigs) - contigs_complete
                    contigs_complete = len(self.completed_contigs)
                    pbar.update(update_number)
            self.write_output(methylation_lines)
            if len(self.completed_contigs) == len(self.contigs) and self.return_queue.empty():
                self.methylation_calling = False
        if self.verbose:
            pbar.close()
        for out in self.output_objects.values():
            out.close()
        self.print_stats()

    def write_output(self, methylation_lines: List[Tuple[Union[int, float, str]]]):
        """Give a list of methylation call dicts, output formatted line

        Params:

        * *methylation_lines (list)*: list of dict containing methylation call information
        """
        # write wig contig designation
        if methylation_lines:
            for meth_tuple in methylation_lines:
                # unpack methylation data
                meth_line = self.unpack_meth_line(meth_tuple)
                # collect methylation stats
                self.collect_stats(meth_line)
                # write ATCGmap line
                if self.ATCGmap:
                    self.write_line(self.output_objects['ATCGmap'], self.format_atcg(meth_line))
                # if methylation level greater than or equal to min_read_depth output CGmap and wig lines
                if meth_line['all_cytosines'] >= self.min_read_depth:
                    self.write_line(self.output_objects['CGmap'], self.format_cgmap(meth_line))

    @staticmethod
    def unpack_meth_line(meth_line: Tuple[Union[int, float, str]]) -> Dict[str, Union[int, float, str]]:
        meth_keys = ['nucleotide', 'meth_cytosines', 'unmeth_cytosines', 'all_cytosines', 'meth_level',
                     'forward_counts', 'reverse_counts', 'pos', 'chrom', 'context', 'subcontext']
        return {key: value for key, value in zip(meth_keys, meth_line)}

    @staticmethod
    def format_atcg(meth_line: Dict[str, Union[int, float, str]]) -> str:
        """ATCGmap line formatting
        """
        ATCGmap_line = f'{meth_line["chrom"]}\t{meth_line["nucleotide"]}\t{meth_line["pos"]}\t{meth_line["context"]}' \
                       f'\t{meth_line["subcontext"]}\t{meth_line["forward_counts"]}\t{meth_line["reverse_counts"]}' \
                       f'\t{meth_line["meth_level"]}\n'
        return ATCGmap_line

    @staticmethod
    def format_cgmap(meth_line: Dict[str, Union[int, float, str]]) -> str:
        """CGmap line formatting
        """
        CGmap_line = f'{meth_line["chrom"]}\t{meth_line["nucleotide"]}\t{meth_line["pos"]}' \
                     f'\t{meth_line["context"]}\t{meth_line["subcontext"]}\t{meth_line["meth_level"]}' \
                     f'\t{meth_line["meth_cytosines"]}\t{meth_line["all_cytosines"]}\n'
        return CGmap_line

    @property
    def get_output_objects(self):
        output_objects = {}
        output_path = self.input_file
        if self.output_prefix:
            output_path = self.output_prefix
        if self.text_output:
            output_objects['CGmap'] = open(f'{output_path}.CGmap', 'w')
            if self.ATCGmap:
                output_objects['ATCGmap'] = open(f'{output_path}.ATCGmap', 'w')
        else:
            output_objects['CGmap'] = io.BufferedWriter(gzip.open(f'{output_path}.CGmap.gz', 'wb'))
            if self.ATCGmap:
                output_objects['ATCGmap'] = io.BufferedWriter(gzip.open(f'{output_path}.ATCGmap.gz', 'wb'))
        return output_objects

    def write_line(self, output_object, line: str):
        """ Outputs line, and encodes if necessary

        Params:

        * *output_object (TextIO/GZipIO)*: output object
        * *line (str)*: formatted line to write
        """
        if self.text_output:
            output_object.write(line)
        else:
            output_object.write(line.encode('utf-8'))

    def collect_stats(self, meth_line):
        """Collect global methylation statistics"""
        if meth_line['subcontext'] == 'CG':
            self.methylation_stats['CG_all'] += meth_line['all_cytosines']
            self.methylation_stats['CG_meth'] += meth_line['meth_cytosines']
        else:
            self.methylation_stats['CH_all'] += meth_line['all_cytosines']
            self.methylation_stats['CH_meth'] += meth_line['meth_cytosines']

    def print_stats(self):
        """Print global methylation statistics"""
        try:
            cpg_percentage = (self.methylation_stats["CG_meth"] / self.methylation_stats["CG_all"]) * 100
        except ZeroDivisionError:
            print('Warning! No SAM reads passed quality control metrics')
            raise MethylationCallingError
        try:
            ch_percentage = (self.methylation_stats["CH_meth"] / self.methylation_stats["CH_all"]) * 100
        except ZeroDivisionError:
            ch_percentage = 0.000
        print(f'Methylated CpG Cytosines: {self.methylation_stats["CG_meth"]}')
        print(f'Total Observed CpG Cytosines: {self.methylation_stats["CG_all"]}')
        print(f'Methylated / Total Observed CpG Cytosines: {cpg_percentage:.2f}%')
        print(f'Methylated CH Cytosines: {self.methylation_stats["CH_meth"]}')
        print(f'Total Observed CH Cytosines: {self.methylation_stats["CH_all"]}')
        print(f'Methylated / Total Observed CH Cytosines: {ch_percentage:.2f}%')
