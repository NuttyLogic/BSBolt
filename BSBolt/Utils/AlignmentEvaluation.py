from collections import defaultdict
from typing import Dict, List, Tuple
import numpy as np
import pysam
from tqdm import tqdm
from BSBolt.Utils.FastqIterator import OpenFastq


def open_alignment_file(alignment_file_path: str = None):
    if alignment_file_path.endswith('.sam'):
        alignment_file = pysam.AlignmentFile(alignment_file_path, 'r')
    else:
        alignment_file = pysam.AlignmentFile(alignment_file_path, 'rb')
    for aligned_segment in alignment_file.fetch(until_eof=True):
        yield aligned_segment


def parse_alignment_comment(read_comment):
    chrom, start, end, cigar, ref_strand = read_comment.replace('+', '').strip().split(':')
    return dict(chrom=chrom, start=int(start),
                end=int(end), cigar=cigar, ref_strand=ref_strand)


def get_read_reference_info(fastq_files: List[str] = None):
    reference_info = {}
    fastq_iterators = (OpenFastq(file) for file in fastq_files)
    for fastq in fastq_iterators:
        for line in fastq:
            read_name = line[0].replace('@', '').strip()
            read_info = parse_alignment_comment(line[2].strip())
            reference_info[read_name] = read_info
    return reference_info


class AlignmentEvaluator:
    """Evaluate alignment against simulated bisulfite sequencing data.

    Params:

    * *duplicated_regions (dict)*: regions duplicated in the simulation reference, [None]
    * *matching_target_prop (float)*: proportion of alignment that most overlap with target region
                                     for a valid alignment to be called, [0.95]
    """

    def __init__(self, duplicated_regions: Dict[str, Tuple[int, int]] = None, matching_target_prop: float = 0.95,
                 verbose: bool = False):
        self.duplicated_regions = duplicated_regions
        self.matching_target_prop = matching_target_prop
        self.verbose = verbose

    def evaluate_alignment(self, alignment_file: str, fastq_files: List[str] = None) -> Dict[str, int]:
        """

        Params:

        * *alignment_file (str)*: path to alignment file
        * *fastq_files (list)*: list of paths to fastq files

        Returns:

        * *alignment_evaluations (dict)*: target alignment stats
        """
        alignment_evaluations = defaultdict(int)
        reads_observed = set()
        reference_info = get_read_reference_info(fastq_files=fastq_files)
        for alignment in tqdm(open_alignment_file(alignment_file), disable=True if not self.verbose else False):
            # read shouldn't have read pair info but added for safety and tools that don't remove
            read_name = alignment.qname.split('/')[0]
            if alignment.is_read1:
                read_name = f'{read_name}/1'
            else:
                read_name = f'{read_name}/2'
            reads_observed.add(read_name)
            alignment_evaluations['ObservedAlignments'] += 1
            alignment_info = reference_info[read_name]
            if not alignment.is_unmapped:
                dup_region = False
                if alignment_info['chrom'] in self.duplicated_regions:
                    duplicate_range = self.duplicated_regions[alignment_info['chrom']]
                    if duplicate_range[0] <= alignment_info['start'] <= duplicate_range[1]:
                        dup_region = True
                chrom_match, matching_prop = self.assess_alignment(alignment, alignment_info)
                on_target = 'On' if chrom_match and matching_prop >= self.matching_target_prop else 'Off'
                secondary = 'Sec' if alignment.is_secondary else 'Prim'
                pair = 'PropPair' if alignment.is_proper_pair else 'Discord'
                dup = 'Dup' if dup_region else 'NoDup'
                alignment_evaluations[f'{on_target}_{secondary}_{pair}_{dup}'] += 1
            else:
                alignment_evaluations['UnalignedAlignments'] += 1
        return alignment_evaluations

    @staticmethod
    def assess_alignment(alignment: pysam.AlignedSegment, alignment_info: Dict):
        """ Compare alignment against simulated reference"""
        chrom_match = alignment.reference_name == alignment_info['chrom']
        # assess reference bases that match between the two reads
        matching_pos = np.array(alignment.get_reference_positions(full_length=False))
        base_range = (matching_pos >= alignment_info['start']) & (matching_pos <= alignment_info['end'])
        matching_prop = sum(base_range) / len(alignment_info['cigar'])
        return chrom_match, matching_prop
