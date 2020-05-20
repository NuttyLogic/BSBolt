import re
from typing import List, Tuple

from BSBolt.Index.RRBSCutSites import ProcessCutSites
from BSBolt.Index.IndexOutput import IndexOutput
from BSBolt.Utils.FastaIterator import OpenFasta
from BSBolt.Utils.UtilityFunctions import get_external_paths


class RRBSBuild:
    """Format reference sequence inputs for processing by BWA. In silico digests reference sequence and
    return mappable regions that are within the fragment boundary. Fragments are relative to the restriction cut site
    if provided, or the complete restriction sequence is considered as part of the mappable fragment.:

    Params:

    * *reference_file (str)*: path to reference file in fasta format
    * *genome_database (str)*: directory to output processed datafiles
    * *block_size (int)*: bwa indexing block size, increases indexing speed but increases memory consumption
    * *lower_bound (int)*: smallest mappable fragment size
    * *upper_bound (int)*: largest mappable fragment size
    * *cut_format (str)*: Comma separated list of restriction sites, - represent cut break
    * *ignore_alt (bool)*: ignore alt contigs when constructing alignment index

    Usage:
    ```python
    index = RRBSBuild(**kwargs)
    index.generate_rrbs_database()
    ```
        """

    def __init__(self, reference_file: str = None, genome_database: str = None,
                 lower_bound: int = 30, upper_bound: int = 500, cut_format: str = 'C-CGG',
                 block_size: int = None, ignore_alt: bool = False):
        bwa_path, _ = get_external_paths()
        self.reference_file = OpenFasta(fasta=reference_file)
        self.index_output = IndexOutput(genome_database=genome_database, block_size=block_size)
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.ignore_alt = ignore_alt
        self.cut_sites = ProcessCutSites(cut_format=cut_format)
        self.mappable_regions = []
        self.contig_size_dict = {}

    def generate_rrbs_database(self):
        """ Wrapper for class functions to process and build mapping indices.
        """
        contig_id = None
        contig_sequence = []
        # return true, line if '>' in fasta line
        for is_contig_label, sequence in self.reference_file:
            if is_contig_label:
                if contig_id:
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
        self.index_output.build_index()
        # output mappable regions in bed format
        self.index_output.output_mappable_regions(self.mappable_regions)
        # output contig size index
        self.index_output.output_contig_sequence('genome_index', self.contig_size_dict)

    def process_contig_region(self, contig_id: str, contig_sequence: List[str]):
        """Given a contig_id will output a pickle file of the whole sequence and output a masked version of the
        the sequence where only mappable regions are reported.

        Params:

        * *contig_id (str)*: contig label
        * *contig_sequence (list)*: a list of of string containing DNA Sequence
        """
        if self.ignore_alt and 'alt' in contig_id.lower():
            return
        # join contig_squence to get ease downstream processing
        contig_str: str = ''.join(contig_sequence)
        # save contig size to dict
        self.contig_size_dict[contig_id] = len(contig_str)
        # serialize contig sequence
        self.index_output.output_contig_sequence(contig_id, contig_str)
        # retrieve list of mappable regions
        mappable_regions: list = self.process_rrbs_sequence(contig_str)
        # check if contig regions empty, if so unmask first 80 bases for indexing
        if not mappable_regions:
            mappable_regions.append((1, 80))
        # get masked contig sequence
        masked_contig_sequence: str = self.mask_contig(contig_str, mappable_regions)
        # perform sanity check, if AssertionError the region designation process is bad
        assert len(contig_str) == len(masked_contig_sequence), 'Contig Length != Masked Contig Length'
        # write contig sequence to output files
        self.index_output.write_contig_sequence(contig_id, masked_contig_sequence)

    def process_rrbs_sequence(self, contig_str: str) -> List[Tuple[int, int]]:
        """Designate mappable regions by finding all occurrences of the restriction site string in the passed DNA
        sequence. Merge restriction map into regions by considering pairs of downstream and upstream restriction sites
        that pass the size limits.

        Params:

        * *contig_str (str)*: STR of continuous DNA sequence

        Returns:

        * *mappable_regions (list)*: List of tuples the contain the start and end position of fragments that
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
                # extend mappable boundaries by two to assess methylation context relative to reference sequence
                mappable_regions.append((down_start - 2, up_end + 2))
        return mappable_regions

    def mask_contig(self, contig_str: str, mappable_regions: List[Tuple[int, int]]) -> str:
        """Given a list of mappable regions, if cut site isn't designated merges mappable fragments, returns a string
                of DNA sequence with masked unmappable regions

                Params:

                * *contig_id (str)*: contig label
                * *contig_str (str)*: str of DNA sequence
                * *mappable_regions (list)*: list of mappable regions

                Returns:

                * *masked_contig_sequence (str)*: str of DNA sequence with un-mappable regions masked
        """
        contig_mappable_regions = iter(mappable_regions)
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
