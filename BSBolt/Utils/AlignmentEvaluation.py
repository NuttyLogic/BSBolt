from collections import defaultdict
from typing import Dict, List, NamedTuple, Tuple
import numpy as np
import pysam
from tqdm import tqdm


def open_alignment_file(alignment_file_path: str = None):
    if alignment_file_path.endswith('.sam'):
        alignment_file = pysam.AlignmentFile(alignment_file_path, 'r')
    else:
        alignment_file = pysam.AlignmentFile(alignment_file_path, 'rb')
    for aligned_segment in alignment_file.fetch(until_eof=True):
        yield aligned_segment


def get_reference_alignments(alignment_file_path):
    alignments = {}
    for alignment in open_alignment_file(alignment_file_path):
        alignment_id = f'{alignment.query_name.split("/")[0]}_1' if alignment.is_read1 else f'{alignment.query_name.split("/")[0]}_2'
        alignments[alignment_id] = alignment.to_string()
    return alignments


class AlignmentEvaluator:

    def __init__(self, duplicated_regions: Dict[str, Tuple[int, int]] = None, matching_target_prop: float = 0.95,
                 verbose: bool = False):
        self.duplicated_regions = duplicated_regions
        self.matching_target_prop = matching_target_prop
        self.verbose = verbose

    def evaluate_alignment(self, alignment_file: str) -> Dict[str, int]:
        alignment_evaluations = defaultdict(int)
        reads_observed = set()
        for alignment in tqdm(open_alignment_file(alignment_file), disable=True if not self.verbose else False):
            read_name = alignment.qname.split('/')[0]
            reads_observed.add(read_name)
            alignment_evaluations['ObservedAlignments'] += 1
            alignment_info = self.parse_alignment_name(read_name)
            if not alignment.is_unmapped:
                dup_region = False
                if alignment_info['chrom'] in self.duplicated_regions:
                    duplicate_range = self.duplicated_regions[alignment_info['chrom']]
                    if duplicate_range[0] <= alignment_info['start'] <= duplicate_range[1]:
                        dup_region = True
                chrom_match, matching_prop = self.assess_alignment(alignment)
                on_target = 'On' if chrom_match and matching_prop >= self.matching_target_prop else 'Off'
                secondary = 'Sec' if alignment.is_secondary else 'Prim'
                pair = 'PropPair' if alignment.is_proper_pair else 'Discord'
                dup = 'Dup' if dup_region else 'NoDup'
                alignment_evaluations[f'{on_target}_{secondary}_{pair}_{dup}'] += 1
            else:
                alignment_evaluations['UnalignedAlignments'] += 1
        return alignment_evaluations

    @staticmethod
    def parse_alignment_name(read_name):
        read_id, chrom, start, end, cigar, ref_strand = read_name.replace('@', '').strip().split(':')
        return dict(read_id=read_id, chrom=chrom, start=int(start), end=int(end), cigar=cigar, ref_strand=ref_strand)

    @staticmethod
    def assess_alignment(alignment: pysam.AlignedSegment, alignment_info: Dict):
        chrom_match = alignment.reference_name == alignment_info['chrom']
        # assess reference bases that match between the two reads
        matching_pos = alignment.get_reference_positions(full_length=False)
        base_range = np.where(alignment_info['start'] <= matching_pos <= alignment_info['end'])
        matching_prop = sum(base_range) / (alignment_info['cigar'])
        return chrom_match, matching_prop
