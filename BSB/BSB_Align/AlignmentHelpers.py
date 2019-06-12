import re
import pysam
from BSB.BSB_Align.LaunchBowtie2Alignment import Bowtie2Align


def launch_bowtie2_stream(bowtie2_stream_kwargs=None, return_dict=None, genome_database_label=None):
    """ Given a mutliprocessing.mangager dict return sam_line"""
    assert isinstance(bowtie2_stream_kwargs, dict)
    assert isinstance(genome_database_label, str)
    for sam_line in Bowtie2Align(**bowtie2_stream_kwargs):
        sam_label = f'{genome_database_label}_{sam_line["QNAME"]}'
        return_dict[sam_label] = sam_line
    return True


def convert_alpha_numeric_cigar(cigar_string):
    """Convert cigar str representation to cigar tuple representation, MATCH = 0, INS = 1, DEL = 2, SOFTCLIP = 4.
    Arguments:
     cigar_string (str): cigar string
    Returns:
         cigar_tuple_list (list [tuples]): list of cigar tuples
    """
    cigar_list = re.findall(r'[^\d_]+|\d+', cigar_string)
    cigar_dict = {'M': 0,  'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    cigar_tuple_list = []
    for cigar_character, cigar_count in zip(cigar_list[1::2], cigar_list[0::2]):
        cigar_tuple_list.append((cigar_dict[cigar_character], int(cigar_count)))
    return cigar_tuple_list


def convert_cigar_tuple(cigar_tuple):
    cigar_dict = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X'}
    cigar_str = ''
    for cigar_type, cigar_length in cigar_tuple:
        cigar_str = f'{cigar_str}{cigar_length}{cigar_dict[cigar_type]}'
    return cigar_str


def get_mapping_length(cigar):
    """Logic to reset position for reverse complement of soft clipped read"""
    # start at next cigar type if beginning of read soft clipped
    read_start = None
    read_end = 0
    mapped_region_length = 0
    reference_consumers = {0, 2, 3, 7, 8}
    query_consumers = {0, 1, 4, 7, 8}
    for cigar_type, cigar_length in cigar:
        if cigar_type in reference_consumers:
            read_start = read_end
            mapped_region_length += cigar_length
        if cigar_type in query_consumers:
            read_end += cigar_length
    return mapped_region_length


def write_bam_line(sam_line=None, output_object=None):
    """Convert sam line to bam line and write out"""
    output_values = list(sam_line.values())[0:-3]
    sam1 = '\t'.join(output_values[:-1])
    sam2 = '\t'.join(output_values[-1])
    formatted_read = f'{sam1}\t{sam2}'
    bam_line = pysam.AlignedSegment.fromstring(formatted_read, output_object.header)
    output_object.write(bam_line)
