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

    def evaluate_alignment(self, reference_alignment_file: str, target_alignment_file: str) -> Dict[str, int]:
        reference_alignments = get_reference_alignments(reference_alignment_file)
        alignment_evaluations = defaultdict(int)
        reads_observed = set()
        for alignment in tqdm(open_alignment_file(target_alignment_file), disable=True if not self.verbose else False):
            read_name = alignment.qname.split('/')[0]
            reads_observed.add(read_name)
            alignment_evaluations['ObservedAlignments'] += 1
            if not alignment.is_unmapped:
                ref_pair = 1 if alignment.is_read1 else 2
                ref_alignment = pysam.AlignedSegment.fromstring(reference_alignments[f'{read_name}_{ref_pair}'],
                                                                alignment.header)
                dup_region = False
                if ref_alignment.reference_name in self.duplicated_regions:
                    duplicate_range = self.duplicated_regions[ref_alignment.reference_name]
                    if duplicate_range[0] <= ref_alignment.reference_start <= duplicate_range[1]:
                        dup_region = True
                chrom_match, cigar_match, matching_prop = self.assess_alignment(alignment, ref_alignment)
                on_target = 'On' if chrom_match and matching_prop >= self.matching_target_prop else 'Off'
                secondary = 'Sec' if alignment.is_secondary else 'Prim'
                pair = 'PropPair' if alignment.is_proper_pair else 'Discord'
                dup = 'Dup' if dup_region else 'NoDup'
                alignment_evaluations[f'{on_target}_{secondary}_{pair}_{dup}'] += 1
            else:
                alignment_evaluations['UnalignedAlignments'] += 1
        alignment_evaluations['UnobservedReads'] += len(reference_alignments) / 2 - len(reads_observed)
        alignment_evaluations['TotalReads'] = len(reference_alignments)
        return alignment_evaluations

    @staticmethod
    def assess_alignment(alignment: pysam.AlignedSegment, ref_alignment: pysam.AlignedSegment):
        chrom_match = alignment.reference_name == ref_alignment.reference_name
        cigar_match = alignment.cigarstring == ref_alignment.cigarstring
        # assess reference bases that match between the two reads
        ref_positions = ref_alignment.get_reference_positions(full_length=False)
        matching_pos = np.isin(alignment.get_reference_positions(full_length=False), ref_positions, assume_unique=True)
        matching_prop = sum(matching_pos) / len(ref_positions)
        return chrom_match, cigar_match, matching_prop
