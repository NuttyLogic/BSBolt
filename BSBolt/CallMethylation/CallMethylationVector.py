import pickle
import re
import pysam
import numpy as np
from BSBolt.Align.AlignmentHelpers import convert_alpha_numeric_cigar


class CallMethylationVector:
    """
     Experimental
     return vector of methylated sites
    """

    def __init__(self, input_file: str = None, genome_database: str = None,
                 contig: str = None, min_base_quality: int = 0, return_queue=None,
                 cg_only: bool = False, start=None, end=None):
        self.input_file = str(input_file)
        self.input_bam = pysam.AlignmentFile(self.input_file, 'rb')
        self.genome_database = str(genome_database)
        if self.genome_database[-1] != '/':
            self.genome_database = f'{self.genome_database}/'
        self.contig = contig
        self.start = start
        self.mate_flags = self.get_mate_flags
        self.end = end
        self.min_base_quality = min_base_quality
        self.chunk_size = 10000
        self.cg_only = cg_only
        self.return_queue = return_queue
        self.counting_dict = {}

    @property
    def get_mate_flags(self):
        # mate flag check for read pairing
        mate_flags = {67: 131, 323: 387, 115: 179, 371: 435,
                      131: 67, 387: 323, 179: 115, 435: 371}
        return mate_flags

    def call_methylation(self):
        try:
            chrom_seq = self.get_reference_sequence(f'{self.genome_database}{self.contig}.pkl')
        except FileNotFoundError:
            self.return_queue.put([])
            print(f'{self.contig} not found in BSBolt DB, Methylation Calls for {self.contig} skipped. Methylation '
                  f'values should be called using the same DB used for alignment.')
            self.return_queue.put([])
        else:
            self.call_contig(chrom_seq)

    def call_contig(self, chrom_seq):
        """Iterates through bam pileup, calling methylation values if the reference nucleotide is a C or G. Pileup reads
        are buffered and accessed as needed.
        """
        # set search pattern for only CG sites or all Cs
        search_pattern = ('C', 'G')
        if self.cg_only:
            search_pattern = ('CG', 'CG')
        # iterate through pileup
        contig_chunk = []
        methylation_vectors = {}
        for aligned_read in self.input_bam.fetch(contig=self.contig, start=self.start, end=self.end):
            if aligned_read.is_unmapped:
                continue
            # get sequence around pileup site
            reference_seq = chrom_seq[aligned_read.reference_start - 1: aligned_read.reference_end + 2].upper()
            c_search_pattern, reference_nuc, strand, offset = search_pattern[0], 'C', 'watson', 0
            if aligned_read.is_reverse:
                c_search_pattern, reference_nuc, strand = search_pattern[1], 'G', 'crick'
                if self.cg_only:
                    offset = 1
            c_search = [match.start() + offset for match in re.finditer(c_search_pattern, reference_seq)]
            # if reference sequence not found proceed to next sequence
            if not c_search:
                continue
            methylation_calls = self.call_vector(aligned_read, c_search, reference_nuc)
            if not methylation_calls[0]:
                continue
            processed_vector = self.process_methylation_vector(aligned_read, methylation_calls,
                                                               strand, methylation_vectors)
            if processed_vector:
                contig_chunk.append(processed_vector)
            if len(contig_chunk) == self.chunk_size:
                self.return_queue.put((self.contig, contig_chunk))
                contig_chunk = []
        # process reads that didn't have a pair with a observed methylation site
        for call in methylation_vectors.values():
            contig_chunk.append((call['read_name'], call['calls'][1][0], call['calls'][1][-1],
                                 np.asarray(call['calls'][0]), np.asarray(call['calls'][1]), call['flag'],
                                 self.mate_flags[call['flag']], call['strand']))
        if contig_chunk:
            self.return_queue.put((self.contig, contig_chunk))

    def call_vector(self, aligned_read, c_search, reference_nuc):
        methylation_calls = [[], []]
        reference_consumers = {0, 2, 3, 7, 8}
        query_consumers = {0, 1, 4, 7, 8}
        positions = iter(c_search)
        # set relative to genomic position so add reference start and one since capturing the first base
        current_pos = next(positions) + aligned_read.reference_start - 1
        reference_pos = aligned_read.reference_start
        query_sequence = aligned_read.query_sequence
        query_qualities = aligned_read.query_qualities
        print(len(query_sequence))
        print(len(query_qualities))
        query_position = 0
        print(aligned_read.cigartuples)
        while current_pos:
            for cigar_type, cigar_count in aligned_read.cigartuples:
                if cigar_type in reference_consumers and cigar_type in query_consumers:
                    for _ in range(cigar_count):
                        while reference_pos > current_pos:
                            current_pos = self.update_position(positions, aligned_read)
                        if reference_pos == current_pos:
                            print(query_position)
                            if query_qualities[query_position] > self.min_base_quality:
                                query_base = query_sequence[query_position]
                                call_made, methylation_value = self.get_methylation_call(reference_nuc, query_base)
                                if call_made:
                                    methylation_calls[0].append(methylation_value)
                                    methylation_calls[1].append(reference_pos)
                            current_pos = self.update_position(positions, aligned_read)
                        reference_pos += 1
                        query_position += 1
                elif cigar_type in reference_consumers and cigar_type not in query_consumers:
                    for _ in range(cigar_count):
                        reference_pos += 1
                elif cigar_type in query_consumers and cigar_type not in reference_consumers:
                    for _ in range(cigar_count):
                        query_position += 1
        return methylation_calls

    def process_methylation_vector(self, aligned_read, methylation_calls, strand, methylation_vectors):
        if aligned_read.is_proper_pair:
            mate_flag = self.mate_flags[aligned_read.flag]
            mate_pair_label = f'{aligned_read.query_name}_{mate_flag}_{aligned_read.next_reference_start}'
            if mate_pair_label in methylation_vectors:
                paired_calls = methylation_vectors.pop(mate_pair_label)['calls']
                paired_calls[0].extend(methylation_calls[0])
                paired_calls[1].extend(methylation_calls[1])
                paired_calls = self.clean_overlap(paired_calls)
                return (aligned_read.query_name, paired_calls[1][0], paired_calls[1][-1],
                        np.array(paired_calls[0]), np.array(paired_calls[1]),
                        aligned_read.flag, mate_flag, strand)
            else:
                vector_label = f'{aligned_read.query_name}_{aligned_read.flag}_{aligned_read.reference_start}'
                methylation_vectors[vector_label] = {'calls': methylation_calls,
                                                     'read_name': aligned_read.query_name,
                                                     'flag': aligned_read.flag, 'strand': strand}
                return None
        else:
            return (aligned_read.query_name, methylation_calls[1][0], methylation_calls[1][-1],
                    np.array(methylation_calls[0]), np.array(methylation_calls[1]), aligned_read.flag, None, strand)

    @staticmethod
    def clean_overlap(methylation_calls):
        cleaned_calls = []
        cleaned_pos = []
        current_pos = -1
        for meth_call, pos in zip(methylation_calls[0], methylation_calls[1]):
            if pos > current_pos:
                cleaned_calls.append(meth_call)
                cleaned_pos.append(pos)
                current_pos = pos
        return [cleaned_calls, cleaned_pos]

    @staticmethod
    def update_position(positions, aligned_read):
        try:
            current_pos = next(positions)
        except StopIteration:
            return None
        else:
            return current_pos + aligned_read.reference_start - 1

    @staticmethod
    def get_methylation_call(reference_nucleotide, base_call):
        """
        Methylation for a C relative to the sense strand of the reference can only be called using watson reads,
        and G with crick reads
        Arguments
            nucleotide (str): reference nucleotide
            base_call (collections.Counter): watson nucleotides are Uppercase and crick nucleotides lowercase
        Returns:
             methylation call dictionary
        """
        # call cytonsines with watson
        if reference_nucleotide == 'C':
            if base_call == 'C':
                return True, 1
            elif base_call == 'T':
                return True, 0
        elif reference_nucleotide == 'G':
            if base_call == 'G':
                return True, 1
            elif base_call == 'A':
                return True, 0
        return False, None

    @staticmethod
    def get_reference_sequence(path):
        """load serialized reference file from path
        """
        with open(path, 'rb') as genome_file:
            return pickle.load(genome_file)
