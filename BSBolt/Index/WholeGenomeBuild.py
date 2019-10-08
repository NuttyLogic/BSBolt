from BSBolt.Utils.FastaIterator import OpenFasta
from BSBolt.Index.IndexOutput import IndexOutput


class WholeGenomeIndexBuild:
    """Class to build whole genome bisulfite bowtie2 index.
    Keyword Arguments
        reference_file (str): path to reference file in fasta format
        genome_database (str): directory to output processed datafiles
        bowtie2_path (str): path to bowtie2 executable, default = bowtie2 if in path
        bowtie2_threads (int): threads for bowtie2 to use
    Attributes
        self.reference_file (OpenFasts): Instance of OpenFasta to parse input reference file
        self.index_output (IndexOutput): Instance of IndexOutput class to handle file output and external bowtie2
            commands
        self.contig_size_dict (dict): dict of contig lengths
    """

    def __init__(self, reference_file=None, genome_database=None, bowtie2_path=None, bowtie2_threads=1,
                 mappable_regions=None):
        assert isinstance(reference_file, str), 'Reference File Path Invalid, Must be a String'
        self.reference_file = OpenFasta(fasta=reference_file)
        assert isinstance(self.reference_file, OpenFasta)
        self.index_output = IndexOutput(**dict(genome_database=genome_database,
                                               bowtie2_path=bowtie2_path,
                                               bowtie2_threads=bowtie2_threads))
        self.mappable_regions = None
        if mappable_regions:
            self.mappable_regions = self.get_mappable_regions(mappable_regions)
        assert isinstance(self.index_output, IndexOutput)
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
        self.index_output.build_bowtie2_index()

    @staticmethod
    def get_mappable_regions(bed_file):
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

    def mask_contig(self, contig_id, contig_str):
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
