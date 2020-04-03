from typing import Dict, List, Tuple

from BSBolt.Index.IndexOutput import IndexOutput
from BSBolt.Utils.FastaIterator import OpenFasta


class WholeGenomeBuild:
    """Class to build whole genome bisulfite bwa index.

    Params:

    * *reference_file (str)*: path to reference file in fasta format
    * *bwa_path (str)*: path to bwa executable, default = bwa-mem if in path

    Usage:
    ```python
    index = WholeGenomeBuild(**kwargs)
    index.generate_bsb_database()
    ```
    """

    def __init__(self, reference_file: str = None, genome_database: str = None,
                 mappable_regions: str = None):
        self.reference_file = OpenFasta(fasta=reference_file)
        self.index_output = IndexOutput(**dict(genome_database=genome_database))
        self.mappable_regions = None
        if mappable_regions:
            self.mappable_regions = self.get_mappable_regions(mappable_regions)
        self.contig_size_dict = {}

    def generate_bsb_database(self):
        """ Wrapper for class functions to process and build mapping indices.
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
        """Get mappable regions if masking with bed file

        Params:

        * *bed_file (str)*: path to bed file of mappable regions

        Returns:

        * *mappable_regions (list)*: sorted list of mappable regions
        """
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
        """Mask contig sequence outside mappable regions

        Params:

        * *contig_id (str)* contig label
        * *contig_str (str)*: contig sequence

        Returns:

        * *contig_str (str)*: masked sequence
        """
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
