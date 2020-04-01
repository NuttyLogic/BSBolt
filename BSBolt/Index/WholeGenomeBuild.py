from typing import Dict, List, Tuple

from BSBolt.Utils.FastaIterator import OpenFasta
from BSBolt.Index.IndexOutput import IndexOutput


class WholeGenomeIndexBuild:
    """Class to build whole genome bisulfite bwa index.
    Keyword Arguments
        reference_file (str): path to reference file in fasta format
        genome_database (str): directory to output processed datafiles
        bwa_path (str): path to bwa executable, default = bwa-mem if in path
    Attributes
        self.reference_file (OpenFasts): Instance of OpenFasta to parse input reference file
        self.index_output (IndexOutput): Instance of IndexOutput class to handle file output and external bwa
            commands
        self.contig_size_dict (dict): dict of contig lengths
    """

    def __init__(self, reference_file: str = None, genome_database: str = None,
                 bwa_path: str = None, mappable_regions: str = None):
        self.reference_file = OpenFasta(fasta=reference_file)
        self.index_output = IndexOutput(**dict(genome_database=genome_database,
                                               bwa_path=bwa_path))
        self.mappable_regions = None
        if mappable_regions:
            self.mappable_regions = self.get_mappable_regions(mappable_regions)
        self.contig_size_dict = {}

    def generate_bsb_database(self):
        """ Wrapper for class functions to process and build mapping indices. Loops through fasta file and processes
        complete contig sequences.
        """
        contig_sequence = []
        contig_id = None
        # return true, line if '>' in fasta line
        for is_contig_label, sequence in self.reference_file:
            if is_contig_label:
                if contig_id:
                    # join contig sequence
                    contig_str = ''.join(contig_sequence)
                    if self.mappable_regions:
                        contig_str = self.mask_contig(contig_id, contig_str)
                    # serialize contig file
                    self.index_output.output_contig_sequence(contig_id, contig_str)
                    # output contig reference
                    self.index_output.write_contig_sequence(contig_id, contig_str)
                    # get contig length
                    self.contig_size_dict[contig_id] = len(contig_str)
                    contig_sequence = []
                # reset contig_id
                contig_id = sequence.replace('>', '').split()[0]
            else:
                contig_sequence.append(sequence)
        # process remaining contig sequence
        contig_str = ''.join(contig_sequence)
        if self.mappable_regions:
            contig_str = self.mask_contig(contig_id, contig_str)
        self.index_output.output_contig_sequence(contig_id, contig_str)
        self.index_output.write_contig_sequence(contig_id, contig_str)
        self.index_output.database_output.close()
        self.contig_size_dict[contig_id] = len(contig_str)
        self.index_output.output_contig_sequence('genome_index', self.contig_size_dict)
        # launch external commands
        self.index_output.build_index()

    @staticmethod
    def get_mappable_regions(bed_file: str) -> Dict[str, List[Tuple[int, int]]]:
        mappable_regions = {}
        with open(bed_file, 'r') as regions:
            for bed_line in regions:
                if bed_line:
                    chrom, start, end = bed_line.replace('\n', '').split('\t')[0:3]
                    start, end = int(start), int(end)
                    if chrom in mappable_regions:
                        mappable_regions[chrom].append((start, end))
                    else:
                        mappable_regions[chrom] = [(start, end)]
        for region_list in mappable_regions.values():
            region_list.sort(key=lambda x: x[0])
        return mappable_regions

    def mask_contig(self, contig_id: str, contig_str: str) -> str:
        if contig_id in self.mappable_regions:
            contig_mappable_regions = iter(self.mappable_regions[contig_id])
            start, end = next(contig_mappable_regions)
            masked_contig = []
            for position, bp in enumerate(contig_str):
                if position > end:
                    try:
                        start, end = next(contig_mappable_regions)
                    except StopIteration:
                        pass
                if start <= position <= end:
                    masked_contig.append(bp)
                else:
                    masked_contig.append('-')
            return ''.join(masked_contig)
        else:
            return '-' * len(contig_str)
