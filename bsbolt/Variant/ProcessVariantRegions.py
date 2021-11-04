import io
import gzip
import multiprocessing
from typing import Dict, List, Tuple, Union
import pysam
from tqdm import tqdm
from bsbolt.Variant.CallVariants import CallRegionVariation
from bsbolt.Utils.CGmapIterator import OpenCGmap


class VariantCallingError(Exception):
    """Error in methylation calling process"""
    pass


def call_variation(completed_contigs, variant_call_kwargs):
    """ Wrapper to initialize methylation calling
    Arguments:
        completed_contigs (multiprocessing.manager.list): List of contigs with completed variant calls
        variant_call_kwargs (dict): dict of argument for CallRegionVariation
    """
    contig_call = CallRegionVariation(**variant_call_kwargs)
    contig_call.call_variation()
    assert isinstance(contig_call, CallRegionVariation)
    completed_contigs.append(variant_call_kwargs['contig'])


class ProcessVarContigs:
    """
    Multi-threaded contig processing wrapper. Passes thread safe queue to workers and outputs values as CGmap file.

    Params:

   * *input_file (str)*: str to input bam/sam file
   * *genome_database (str)*: str to genome directory
   * *output_prefix (str)*: output prefix for  CGmap
   * *ignore_overlap (bool)*:  ignore overlapping reads, [True]
   * *text_output (bool)*: output compressed or plain text files, [False]
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
                 ignore_overlap: bool = True, text_output: bool = False,
                 min_read_depth: int = 10, max_read_depth: int = 8000, threads: int = 1, verbose: bool = True,
                 min_base_quality: int = 10, min_mapping_quality: int = 10, ignore_orphans: bool = False,
                 bed_output: bool = True, vcf_output: bool = False, output_reference_calls: bool = False,
                 call_region_bed: str = None, min_pval: float = 0.001):
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
        self.call_region_bed = call_region_bed
        self.output_prefix = output_prefix
        self.threads = threads
        self.call_kwargs = dict(input_file=input_file,
                                genome_database=genome_database,
                                ignore_overlap=ignore_overlap,
                                ignore_orphans=ignore_orphans,
                                max_read_depth=max_read_depth,
                                min_base_quality=min_base_quality,
                                min_mapping_quality=min_mapping_quality,
                                min_read_depth=min_read_depth)
        self.min_read_depth = min_read_depth
        self.min_pval = min_pval
        self.bed_output = bed_output
        self.vcf_output = vcf_output
        self.output_ref = output_reference_calls
        self.calling = True
        self.contigs = self.get_contigs
        self.completed_contigs = None
        self.return_queue = None
        self.pool = None
        self.verbose = verbose
        self.output_objects = self.get_output_objects
        self.call_stats = {'Heterozygous': 0, '': 0, 'Homozygous': 0, 'Alt Alleles': 0}

    @property
    def get_contigs(self) -> List[str]:
        """
        get list of contigs in input file, threads set across contigs
        """
        contigs = []
        if self.call_region_bed:
            for line in OpenCGmap(self.call_region_bed):
                contigs.append((line[0], int(line[1]), int(line[2])))
        else:
            contigs = [[contig[0]] for contig in self.input_bam.get_index_statistics() if contig[1] > 0]
        if not contigs:
            print('No reads are mapped\nPlease check alignment file\n')
            raise VariantCallingError
        return contigs

    def process_contigs(self):
        """Launches a processing pool to call methylation values across the input file contigs
        """
        # initialize manager
        manager = multiprocessing.Manager()
        # get return dictionary
        self.return_queue = manager.Queue(maxsize=20)
        self.completed_contigs = manager.list()
        # threads for variant calling, if one thread use thread for calling and watching
        pool_threads = self.threads - 1 if self.threads != 1 else 1
        # start pool
        self.pool = multiprocessing.Pool(processes=pool_threads)
        # for contig call methylation and return values to dict
        for contig in self.contigs:
            contig_kwargs = dict(self.call_kwargs)
            contig_kwargs.update(dict(contig=contig[0], return_queue=self.return_queue))
            if len(contig) > 1:
                contig_kwargs.update(dict(start=contig[1], stop=contig[2]))
            self.pool.apply_async(call_variation,
                                  args=[self.completed_contigs, contig_kwargs],
                                  error_callback=self.variant_calling_error)
        self.pool.close()

    def variant_calling_error(self, error):
        """Raise if exception thrown in methylation calling process"""
        self.variant_calling_error = False
        raise VariantCallingError(error)

    def watch_pool(self):
        """Watch self.return_dict and process return methylation values. Contigs are processed in order so buffer can
        become large if first contig is large, ie. Human Chr1
        """
        contigs_complete = 0
        pbar = None
        if self.verbose:
            pbar = tqdm(total=len(self.contigs), desc='Processing Contigs')
        while self.calling:
            variant_calls: list = self.return_queue.get(block=True)
            if self.verbose:
                if len(self.completed_contigs) != contigs_complete:
                    update_number = len(self.completed_contigs) - contigs_complete
                    contigs_complete = len(self.completed_contigs)
                    pbar.update(update_number)
            self.write_output(variant_calls)
            if len(self.completed_contigs) == len(self.contigs) and self.return_queue.empty():
                self.calling = False
        if self.verbose:
            pbar.close()
        for out in self.output_objects.values():
            out.close()
        self.print_stats()

    def write_output(self, variant_calls: List):
        """Give a list of methylation call dicts, output formatted line

        Params:

        * *methylation_lines (list)*: list of dict containing methylation call information
        """
        # write wig contig designation
        if variant_calls:
            for call in variant_calls:
                # unpack call data
                var_call = self.unpack_var_call(call)
                # ToDo: collect variant call stats 
                if var_call['call_p'] > self.min_pval:
                    continue
                if not self.output_ref:
                    if var_call['genotype'] == var_call['ref_base']:
                        continue
                if self.bed_output:
                    self.write_line(self.output_objects['bed'], self.format_bed(var_call))
                if self.vcf_output:
                    self.write_line(self.output_objects['vcf'], self.format_vcf(var_call))

    @staticmethod
    def unpack_var_call(variant_call):
        var_keys = ('contig', 'pos', 'call_prob', 'call_p', 'call_score', 'ref_base', 'genotype')
        var_call = {key: val for key, val in zip(var_keys, variant_call)}
        if var_call['genotype'][0] == var_call['genotype'][1]:
            var_call['genotype'] = var_call['genotype'][0]
        else:
            var_call['genotype'] = ','.join(list(var_call['genotype']))
        return var_call

    @staticmethod
    def format_bed(variant_call) -> str:
        """bed line formatting
        """
        bed_line = f'{variant_call["contig"]}\t{variant_call["pos"]-1}\t{variant_call["pos"]}\t{variant_call["call_p"]}'\
                   f'\t{variant_call["call_score"]}\t{variant_call["ref_base"]}\t{variant_call["genotype"]}\n'
        return bed_line

    @staticmethod
    def format_vcf(variant_call) -> str:
        """VCF line formatting
        """
        return ' '

    @property
    def get_output_objects(self):
        output_objects = {}
        output_path = self.input_file
        if self.output_prefix:
            output_path = self.output_prefix
        if self.bed_output:
            if self.text_output:
                output_objects['bed'] = open(f'{output_path}.bed', 'w')
            else:
                output_objects['bed'] = io.BufferedWriter(gzip.open(f'{output_path}.bed.gz', 'wb'))
        if self.vcf_output:
            if self.text_output:
                output_objects['vcf'] = open(f'{output_path}.vcf', 'w')
            else:
                output_objects['vcf'] = io.BufferedWriter(gzip.open(f'{output_path}.vcf.gz', 'wb'))
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

    def collect_stats(self, variant_call):
        """Collect global methylation statistics"""
        pass


    def print_stats(self):
        """Print global methylation statistics"""
        pass