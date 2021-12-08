from collections import Counter
import pickle
from typing import Dict, Tuple, Union
import numpy as np
import pysam
from bsbolt.Variant.VariantCall import CallVariant


class CallRegionVariation:
    """
    Params:

    * *input_file (str)*: str to input bam/sam file
    * *genome_database (str)*: str to genome directory
    * *ignore_overlap (bool)*:  ignore overlapping signal, calls methylation using the higher quality base for
                                reads with overlapping alignments, [True]
    * *ignore_oprhans (bool)*: ignore orphaned reads (not properly paired), [True]
    * *max_read_depth (int)*: maximum read depth for pileup, [8000]
    * *contig (str)*: contig to process
    * *start (int)*: contig start position
    * *stop (int)*: contig stop position
    * *min_base_quality (int)*: minimum quality for a base to be reported for methylation calling, [10]
    * *return_queue (Queue.queue)*: results are added to a queue in batches for multi-threaded access
    * *cg_only (bool)*: only return CG sites to queue, [False]
    * *min_mapping_quality (int)*: minimum mapping quality for an alignment to be used for methylation calling, [10]

    Attributes:

    * *self.chunk_size (int)*: Number of sites in chunk added to queue, [10000]

    """

    def __init__(self, input_file: str = None, genome_database: str = None,
                 ignore_overlap: bool = True, ignore_orphans: bool = True, min_read_depth: int = 10,
                 max_read_depth: int = 8000, contig: str = None, min_base_quality: int = 10, return_queue=None,
                 min_mapping_quality: int = 10, start: int = None, stop: int = None):
        self.input_file = str(input_file)
        self.input_bam = pysam.AlignmentFile(self.input_file, mode='rb',
                                             require_index=True)
        self.genome_database = str(genome_database)
        if self.genome_database[-1] != '/':
            self.genome_database = f'{self.genome_database}/'
        self.ignore_overlap = ignore_overlap
        self.ignore_orphans = ignore_orphans
        self.max_read_depth = max_read_depth
        self.min_read_depth = min_read_depth
        self.contig = contig
        self.min_base_quality = min_base_quality
        self.chunk_size = 10000
        self.min_mapping_quality = min_mapping_quality
        self.start = start - 1 if start else None
        self.stop = stop + 1 if start else None
        self.return_queue = return_queue
        self.counting_dict = {}

    def call_variation(self):
        """Run methylation call for contig"""
        try:
            chrom_seq = self.get_reference_sequence(f'{self.genome_database}{self.contig}.pkl')
        except FileNotFoundError:
            print(f'{self.contig} not found in BSBolt DB, Variant Calls for {self.contig} skipped. Methylation '
                  f'values should be called using the same DB used for alignment.')
            self.return_queue.put([])
        else:
            self.call_contig(chrom_seq)

    def call_contig(self, chrom_seq: str):
        """Iterates through bam pileup, calling methylation values if the reference nucleotide is a C or G. Pileup reads
        are buffered and accessed as needed. Process sites are appended to list and returned in chunks.

        Reads flagged as:

        * pcr duplicates (1024), read unmapped (4), read fails platform/vendor quality checks (512) are ignored
        """
        # iterate through pileup
        line_count = 0
        contig_chunk = []
        # flag_require (0) mapped read
        # flag_filter 1540 = pcr duplicates (1024) + read unmapped (4) + read fails platform/vendor quality checks (512)
        caller = CallVariant()
        for pileup_col in self.input_bam.pileup(max_depth=self.max_read_depth,
                                                contig=self.contig,
                                                start=self.start,
                                                stop=self.stop,
                                                ignore_overlaps=self.ignore_overlap,
                                                min_base_quality=self.min_base_quality,
                                                ignore_orphans=self.ignore_orphans,
                                                min_mapping_quality=self.min_mapping_quality,
                                                flag_require=0,
                                                flag_filter=1540):
            if self.start:
                if not self.start < pileup_col.reference_pos < self.stop:
                    continue
            # get sequence around pileup site
            reference_seq = chrom_seq[(pileup_col.reference_pos - 3):(pileup_col.reference_pos + 4)].upper()
            # get nucleotide context
            fivemer = reference_seq[1:-1]
            if len(fivemer) == 5:
                nucleotide = fivemer[2]
                #cg_site = self.check_cg(nucleotide, fivemer)
                # count the pileup read bases, Uppercase watson strand, lowercase crick strand
                try:
                    base_counts = Counter(pileup_col.get_query_sequences(mark_matches=False,
                                                                         mark_ends=False,
                                                                         add_indels=False))
                except AssertionError:
                    # pysam may throw an error if the number of reads in pileup column is above the max
                    continue
                base_counts = np.array([base_counts.get(key, 0) for key in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']])
                if sum(base_counts) < self.min_read_depth:
                    continue
                call = caller.call_variant(base_counts)
                variant_call = dict(chrom=self.contig, pos=pileup_col.reference_pos,
                                    call_prob=call[0], call_p=call[1],
                                    call_score=call[2], ref_base=nucleotide, genotype=call[3])
                contig_chunk.append(tuple(variant_call.values()))
                line_count += 1
                # return chunk for output
                if line_count == self.chunk_size - 1:
                    self.return_queue.put(contig_chunk, block=True)
                    contig_chunk = []
                    line_count = 0
        self.return_queue.put(contig_chunk, block=True)

    @staticmethod
    def check_cg(nucleotide: str, fivemer: str) -> Tuple[str, str]:
        """
        Params:

        * *nucleotide (str)*: 1 nucleotide
        * *fivemer (str)*: 5 nucleotides

        Returns:

        * *cg (bool)*: True if CG site 
        """
        if nucleotide == 'C' and fivemer[2:4] == 'CG':
            return True
        elif nucleotide == 'G' and fivemer[1:3] == 'CG':
            return True
        return False

    @staticmethod
    def get_reference_sequence(path: str) -> str:
        """load serialized reference file from path
        """
        with open(path, 'rb') as genome_file:
            return pickle.load(genome_file)
