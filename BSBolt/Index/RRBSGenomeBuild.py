import re

from BSBolt.Utils.FastaIterator import OpenFasta
from BSBolt.Index.ProcessCutSites import ProcessCutSites
from BSBolt.Index.IndexOutput import IndexOutput


class RRBSGenomeIndexBuild:
    """Class to format reference sequence inputs for processing by bowtie2. In silico digests reference sequence and
    return mappable regions that are within the fragment boundary. Fragments are relative to the restriction cut site
    if provided, or the complete restriction sequence is considered as part of the mappable fragment.
        Keyword Arguments:
            reference_file (str): path to reference file in fasta format
            genome_database (str): directory to output processed datafiles
            bowtie2_path (str): path to bowtie2 executable, default = bowtie2 if in path
            bowtie2_threads (int): threads for bowtie2 to use
            lower_bound (int): smallest mappable fragment size
            upper_bound (int): largest mappable fragment size
            cut_format (str): Comma separated list of restriction sites, - represent cut break
        Attributes:
            self.reference_file (OpenFasta): Instance of OpenFasta to parse input reference file
            self.index_output (IndexOutput): Instance of IndexOutput class to handle file output and external bowtie2
                commands
            self.lower_bound (int): == lower_bound kwarg
            self.higher_bound (int): == high_bound kwarg
            self.cut_sites (ProcessCutSites): ProcessCutSites instance containing recognition sequences and cut offsets
            self.mappable_regions (list): list of regions that are mappable
            self.contig_size_dict (dict): list of contig sizes for downstream use
        """

    def __init__(self, reference_file=None, genome_database=None, bowtie2_path='bowtie2',
                 bowtie2_threads=1, lower_bound=30, upper_bound=500, cut_format='C-CGG'):
        assert isinstance(reference_file, str), 'Reference File Path Invalid, Must be a String'
        assert isinstance(lower_bound, int), 'Bound Invalid, Must be Integer'
        assert isinstance(upper_bound, int), 'Bound Invalid, Must be Integer'
        assert isinstance(cut_format, str), 'Cut Format Invalid, Must be String'
        self.reference_file = OpenFasta(fasta=reference_file)
        assert isinstance(self.reference_file, OpenFasta)
        self.index_output = IndexOutput(**dict(genome_database=genome_database,
                                               bowtie2_path=bowtie2_path,
                                               bowtie2_threads=bowtie2_threads))
        assert isinstance(self.index_output, IndexOutput)
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.cut_sites = ProcessCutSites(cut_format=cut_format)
        self.mappable_regions = []
        self.contig_size_dict = {}

    def generate_rrbs_database(self):
        """ Wrapper for class functions to process and build mapping indices. Loops through fasta file and processes
        complete contig sequences.
        """
        contig_id = None
        contig_sequence = []
        # return true, line if '>' in fasta line
        for is_contig_label, sequence in self.reference_file:
            if is_contig_label:
                if contig_id:
                    # process contig sequence
                    self.process_contig_region(contig_id, contig_sequence)
                    contig_sequence = []
                # set contig_id
                contig_id = sequence.replace('>', '').split()[0]
            else:
                contig_sequence.append(sequence.upper())
        # process remaining sequence after iterator exhausted
        self.process_contig_region(contig_id, contig_sequence)
        self.index_output.database_output.close()
        # launch build commands
        self.index_output.build_bowtie2_index()
        # output mappable regions in bed format
        self.index_output.output_mappable_regions(self.mappable_regions)
        # output contig size index
        self.index_output.output_contig_sequence('genome_index', self.contig_size_dict)

    def process_contig_region(self, contig_id, contig_sequence):
        """Given a contig_id will output a pickle file of the whole sequence and output a masked version of the
        the sequence where only mappable regions are reported.
            Arguments:
                contig_id (str): contig label
                contig_sequence (list): a list of of string containing DNA Sequence
        """
        # join contig_squence to get ease downstream processing
        contig_str: str = ''.join(contig_sequence)
        # save contig size to dict
        self.contig_size_dict[contig_id] = len(contig_str)
        # serialize contig sequence
        self.index_output.output_contig_sequence(contig_id, contig_str)
        # retrieve list of mappable regions
        contig_regions: list = self.process_rrbs_sequence(contig_str)
        # check if contig regions empty, if so unmask first 80 bases to makes sure bt2 indexing is successful
        if not contig_regions:
            contig_regions.append((1, 80))
        # get masked contig sequence
        masked_contig_sequence: str = self.mask_unmappable_sites(contig_id, contig_str, contig_regions)
        # perform sanity check, if AssertionError the region designation process is bad
        assert len(contig_str) == len(masked_contig_sequence), 'Contig Length != Masked Contig Length'
        # write contig sequence to output files
        self.index_output.write_contig_sequence(contig_id, masked_contig_sequence)

    def process_rrbs_sequence(self, contig_str):
        """Designate mappable regions by finding all occurrences of the restriction site string in the passed DNA
        sequence. Merge restriction map into regions by considering pairs of downstream and upstream restriction sites
        that pass the size limits.
            Arguments:
                contig_str (str): STR of continuous DNA sequence
            Returns:
                mappable_regions (list): List of tuples the contain the start and end position of fragments that
                    are with the size limits
        """
        restriction_site_locations = []
        # get the position of the all occurrences of the restriction site pattern in the DNA sequence and add to list
        for restrication_site, offset in self.cut_sites.restriction_site_dict.items():
            restriction_site_locations.extend([(m.start(), offset, restrication_site) for m in
                                               re.finditer(restrication_site, contig_str)])
        # sort list so fragments are in order or start position
        restriction_site_locations.sort(key=lambda x: x[0])
        mappable_regions = []
        # iterate through list to get the cut site plus upstream cut site, will terminate at the second to last site
        for downstream_site, upstream_site in zip(restriction_site_locations, restriction_site_locations[1:]):
            # set default downstream site that includes the restriction site
            down_start = downstream_site[0]
            # if cut site indicated shift start position
            if downstream_site[1]:
                down_start = downstream_site[0] + downstream_site[1]
            # set site end to include restriction site
            up_end = upstream_site[0] + len(upstream_site[2])
            # correct fragment if cut site provided
            if upstream_site[1]:
                up_end = upstream_site[0] + upstream_site[1]
            site_range = up_end - down_start
            # check fragment length
            if self.lower_bound <= site_range <= self.upper_bound:
                mappable_regions.append((down_start, up_end))
        return mappable_regions

    def mask_unmappable_sites(self, contig_id, contig_str, contig_regions, masking_nucleotide='-'):
        """Given a list of mappable regions, if cut site isn't designated merges mappable fragments, returns a string
        of DNA sequence with masked unmappable regions
        Arguments:
            contig_id (str): contig label
            contig_str (str): str of DNA sequence
            contig_regions (list): list of mappable regions
            masking_nucleotide (str): str to use for nucleotide maksing
        Keyword Arguments:
            masking_nucleotide (str): default = '-', character used to masked DNA sequence
        Returns:
            masked_contig_sequence (str): str of DNA sequence with un-mappable regions masked
        """
        masked_contig_sequence = []
        # initialize last index at 0
        ending_index = 0
        for region in contig_regions:
            # if the regions aren't overlapping used region start site, else use the ending_index
            masked_sequence = (region[0] - ending_index) * masking_nucleotide
            mappable_region = contig_str[region[0]:region[1]]
            masked_contig_sequence.append(masked_sequence)
            masked_contig_sequence.append(mappable_region)
            ending_index = region[1]
            # add mapping region info to self.mappable_regions list in bed format
            self.mappable_regions.append(f'{contig_id}\t{region[0]}\t{region[1] - 1}\t{mappable_region}\n')
        # append remaining contig sequence as masked, by definition the last part can't be a mappable region
        masked_contig_sequence.append(len(contig_str[ending_index:]) * masking_nucleotide)
        return ''.join(masked_contig_sequence)
