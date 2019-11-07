import re
from typing import List, Tuple
import pysam
from BSBolt.Align.LaunchBowtie2Alignment import Bowtie2Align


def convert_alpha_numeric_cigar(cigar_string):
    """Convert cigar str representation to cigar tuple representation, MATCH = 0, INS = 1, DEL = 2, SOFTCLIP = 4.
    Arguments:
     cigar_string (str): cigar string
    Returns:
         cigar_tuple_list (list [tuples]): list of cigar tuples
    """
    cigar_list = re.findall(r'[^\d_]+|\d+', cigar_string)
    cigar_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    cigar_tuple: List[Tuple[int, int]] = []
    for cigar_character, cigar_count in zip(cigar_list[1::2], cigar_list[0::2]):
        cigar_tuple.append((cigar_dict[cigar_character], int(cigar_count)))
    return cigar_tuple


def convert_cigar_tuple(cigar_tuple):
    """Convert cigar tuple to cigar string, assumes cigar tuple is from strand being reverse complemented"""
    cigar_dict = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X'}
    cigar_str = ''
    for cigar_type, cigar_length in cigar_tuple:
        cigar_str = f'{cigar_str}{cigar_length}{cigar_dict[cigar_type]}'
    return cigar_str


def get_mapping_length(cigar):
    """Logic to reset position for reverse complement of soft clipped read"""
    # start at next cigar type if beginning of read soft clipped
    read_end = 0
    mapped_region_length = 0
    reference_consumers = {0, 2, 3, 7, 8}
    query_consumers = {0, 1, 4, 7, 8}
    for cigar_type, cigar_length in cigar:
        if cigar_type in reference_consumers:
            mapped_region_length += cigar_length
        if cigar_type in query_consumers:
            read_end += cigar_length
    return mapped_region_length


def write_bam_line(sam_line=None, output_object=None):
    """Convert sam line to bam line and write out"""
    sam_line_1 = f'{sam_line["QNAME"]}\t{sam_line["FLAG"]}\t{sam_line["RNAME"]}\t{sam_line["POS"]}\t'
    sam_line_2 = f'{sam_line["MAPQ"]}\t{sam_line["CIGAR"]}\t{sam_line["RNEXT"]}\t{sam_line["PNEXT"]}\t'
    sam_line_3 = f'{sam_line["TLEN"]}\t{sam_line["SEQ"]}\t{sam_line["QUAL"]}\t'
    sam_tags = '\t'.join(sam_line['SAM_TAGS'])
    formatted_read = f'{sam_line_1}{sam_line_2}{sam_line_3}{sam_tags}'
    bam_line = pysam.AlignedSegment.fromstring(formatted_read, output_object.header)
    output_object.write(bam_line)


def bowtie2_alignment(bowtie2_kwargs, read_pipe, paired_end=False, chunk_size=10000):
    alignment = iter(Bowtie2Align(bowtie2_kwargs))
    sam_reads, read_name = [], None
    read_chunk, read_count = [], 0
    while True:
        try:
            sam_read = next(alignment)
        except StopIteration:
            break
        if not read_name:
            read_name = sam_read['QNAME']
        elif sam_read['QNAME'] != read_name:
            read_chunk.append(sam_reads)
            sam_reads, read_name = [], sam_read['QNAME']
        if paired_end:
            sam_read_paired = next(alignment)
            assert sam_read['QNAME'] == sam_read_paired['QNAME'], '.fastq files are not paired, sort .fastq files'
            sam_reads.append((sam_read, sam_read_paired))
        else:
            sam_reads.append(tuple(sam_read))
        read_count += 1
        if read_count == chunk_size:
            read_pipe.send(read_chunk)
            read_chunk, read_count = [], 0
    read_pipe.send(read_chunk)
