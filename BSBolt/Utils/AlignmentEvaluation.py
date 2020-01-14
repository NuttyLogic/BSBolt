from collections import namedtuple
from typing import Dict, List, NamedTuple, Tuple
import numpy as np
import pysam

alignment_type = NamedTuple('alignment', [('name', str), ('flag', str), ('ref_name', str), ('ref_pos', str),
                                          ('map_quality', str), ('cigar', str), ('next_ref_name', str),
                                          ('next_ref_pos', str), ('length', str), ('seq', str), ('qual', str),
                                          ('tags', List[str]), ('alignment_positions', np.ndarray)])


def open_alignment_file(alignment_file_path: str = None):
    if alignment_file_path.endswith('.sam'):
        alignment_file = pysam.AlignmentFile(alignment_file_path, 'r')
    else:
        alignment_file = pysam.AlignmentFile(alignment_file_path, 'rb')
    alignment = namedtuple('alignment', ['name', 'flag', 'ref_name', 'ref_pos', 'map_quality', 'cigar',
                                         'next_ref_name', 'next_ref_pos', 'length', 'seq', 'qual', 'tags',
                                         'alignment_positions'])
    for line in alignment_file.fetch(until_eof=True):
        mapping_positions = line.get_reference_positions()
        read_info = list(line.to_dict().values()) + [mapping_positions]
        yield alignment(*read_info)


def get_alignments(alignment_file_path):
    alignments = {}
    for alignment in open_alignment_file(alignment_file_path):
        if alignment.name in alignments:
            alignments[alignment.name].append(alignment)
        else:
            alignments[alignment.name] = [alignment]
    return alignments


class AlignmentEvaluator:

    def __init__(self, duplicated_regions: Dict[str, Tuple[int, int]] = None, matching_target_prop: float = 0.95):
        self.first_pair_flags = {'65', '67', '69', '73', '77', '81', '83', '97', '99', '101',
                                 '107', '113', '115', '321', '323', '325', '329', '333', '337',
                                 '339', '353', '355', '357', '363', '369', '371'}
        self.discordant_flags = ['65', '81', '113', '129', '145', '161', '177', '73', '89', '137', '153']
        self.discordant_flags = self.discordant_flags + [str(int(flag) + 256) for flag in self.discordant_flags]
        self.proper_pair_flags = ['67', '83', '99', '115', '131', '147', '163', '179']
        self.proper_pair_flags = self.proper_pair_flags + [str(int(flag) + 256) for flag in self.proper_pair_flags]
        self.mapping_flags = set(self.discordant_flags + self.proper_pair_flags)
        self.mixed_unmapped = {'69', '101', '133', '165'}
        self.duplicated_regions = duplicated_regions
        self.matching_target_prop = matching_target_prop

    def evaluate_alignment(self, reference_alignments, target_alignments):
        alignment_evaluations = self.asses_alignments(reference_alignments, target_alignments)
        return self.summarize_evaluations(alignment_evaluations)

    @staticmethod
    def summarize_evaluations(alignment_evaluations: List[Dict]):
        summary = dict(total_reads=0, on_target_paired=0, on_target_single=0, off_target_paired=0,
                       off_target_single=0, supplementary_dup=0, supplementary=0, unaligned=0, total_alignments=0)
        for alignment_group in alignment_evaluations:
            # on target alignments are alignments without supplementary alignments or occur within a duplicated region
            summary['total_reads'] += 1
            summary['on_target_paired'] += alignment_group['on_target_paired']
            summary['on_target_single'] += alignment_group['on_target_single']
            summary['unaligned'] += alignment_group['unaligned']
            if not alignment_group['dup_region']:
                summary['off_target_paired'] += alignment_group['off_target_paired']
                summary['off_target_single'] += alignment_group['off_target_single']
                summary['supplementary'] += alignment_group['observed'] - 2
                summary['total_alignments'] += alignment_group['observed']
            else:
                # account for duplicated region
                off_target_paired = alignment_group['off_target_paired'] - 2
                summary['off_target_paired'] += off_target_paired if off_target_paired > -1 else 0
                off_target_single = alignment_group['off_target_single'] - 2
                summary['off_target_single'] += off_target_single if off_target_single > -1 else 0
                summary['supplementary_dup'] += alignment_group['observed'] - 2
                summary['total_alignments'] += alignment_group['observed'] - 2
        return summary

    def asses_alignments(self,
                         reference_alignments: Dict[str, List[alignment_type]],
                         target_alignments: Dict[str, List[alignment_type]]):
        alignment_evaluations = []
        for read_name, read_alignments in target_alignments.items():
            ref_alignments = self.format_reference(reference_alignments[read_name])
            dup_region = False
            for alignment in ref_alignments.values():
                if alignment.ref_name in self.duplicated_regions:
                    duplicate_range = self.duplicated_regions[alignment.ref_name]
                    if duplicate_range[0] <= int(alignment.ref_pos) <= duplicate_range[1]:
                        dup_region = True
                        break
            alignment_quality = dict(on_target_paired=0, off_target_paired=0, on_target_single=0, off_target_single=0,
                                     unaligned=0, unaligned_mixed=0, observed=0, dup_region=dup_region)
            for alignment in read_alignments:
                alignment_quality['observed'] += 1
                ref_pair = 1 if alignment.flag in self.first_pair_flags else 2
                chrom_match, cigar_match, matching_prop = self.assess_alignment(alignment, ref_alignments[ref_pair])
                if alignment.flag in self.mapping_flags:
                    paired = alignment.flag in self.proper_pair_flags
                    on_target = chrom_match and matching_prop >= self.matching_target_prop
                    if on_target and paired:
                        alignment_quality['on_target_paired'] += 1
                    elif on_target and not paired:
                        alignment_quality['on_target_single'] += 1
                    elif not on_target and paired:
                        alignment_quality['off_target_paired'] += 1
                    else:
                        alignment_quality['off_target_single'] += 1
                else:
                    if alignment.flag in self.mixed_unmapped:
                        alignment_quality['unaligned_mixed'] += 1
                    else:
                        alignment_quality['unaligned'] += 1
            alignment_evaluations.append(alignment_quality)
        return alignment_evaluations

    def format_reference(self, ref_alignments: List[alignment_type]):
        formatted_ref = {1: None, 2: None}
        for alignment in ref_alignments:
            if alignment.flag in self.first_pair_flags:
                formatted_ref[1] = alignment
            else:
                formatted_ref[2] = alignment
        return formatted_ref

    @staticmethod
    def assess_alignment(alignment: alignment_type, ref_alignment: alignment_type):
        chrom_match = alignment.ref_name == ref_alignment.ref_name
        cigar_match = alignment.cigar == ref_alignment.cigar
        # assess reference bases that match between the two reads
        matching_pos = np.isin(alignment.alignment_positions, ref_alignment.alignment_positions, assume_unique=True)
        matching_prop = sum(matching_pos) / len(ref_alignment.alignment_positions)
        return chrom_match, cigar_match, matching_prop
