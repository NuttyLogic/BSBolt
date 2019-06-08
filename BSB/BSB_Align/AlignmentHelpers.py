import re
import pysam
from BSB.BSB_Align.LaunchBowtie2Alignment import Bowtie2Align
from BSB.BSB_Align.StreamTabFormat import StreamTab


def launch_bowtie2_stream(bowtie2_stream_kwargs=None, return_dict=None, genome_database_label=None):
    """ Given a mutliprocessing.mangager dict return sam_line"""
    assert isinstance(bowtie2_stream_kwargs, dict)
    assert isinstance(genome_database_label, str)
    for sam_line in Bowtie2Align(**bowtie2_stream_kwargs):
        sam_label = f'{genome_database_label}_{sam_line["QNAME"]}'
        return_dict[sam_label] = sam_line
    return True


def launch_unaltered_tab_stream(tab_kwargs=None, return_dict=None, read_list=None):
    """ Stream tab format without altering base"""
    for tab_dict in StreamTab(**tab_kwargs):
        return_dict.update(tab_dict)
        read_list.extend(list(tab_dict.keys()))


def convert_alpha_numeric_cigar(cigar_string):
    """Convert cigar str representation to cigar tuple representation, MATCH = 0, INS = 1, DEL = 2, SOFTCLIP = 4.
    Arguments:
     cigar_string (str): cigar string
    Returns:
         cigar_tuple_list (list [tuples]): list of cigar tuples
    """
    cigar_list = re.findall(r'[^\d_]+|\d+', cigar_string)
    cigar_dict = {'M': 0,  'I': 1, 'D': 2, 'S': 4, 'N': 0}
    cigar_tuple_list = []
    for cigar_character, cigar_count in zip(cigar_list[1::2], cigar_list[0::2]):
        cigar_tuple_list.append((cigar_dict[cigar_character], int(cigar_count)))
    return cigar_tuple_list


def get_length_stats(cigar):
    """Logic to reset position for reverse complement of soft clipped read"""
    # start at next cigar type if beginning of read soft clipped
    read_start = cigar[0][1] if cigar[0][0] == 4 else 0
    read_end = int(read_start)
    mapped_region_length = 0
    for cigar_type, cigar_length in cigar:
        if cigar_type == 0:
            read_end += cigar_length
            mapped_region_length += cigar_length
        elif cigar_type == 1:
            read_end += cigar_length
        elif cigar_type == 2:
            mapped_region_length += cigar_length
    return read_start, read_end, mapped_region_length


def launch_bowtie2_mapping(bowtie2_stream_kwargs=None, output_path=None):
    """Write output of Bowtie2Align instance"""
    assert isinstance(bowtie2_stream_kwargs, dict)
    output_object = open(f'{output_path}.sam.temp', 'w')
    for sam_line in Bowtie2Align(**bowtie2_stream_kwargs):
        write_sam_line(sam_line=sam_line, output_object=output_object)
    output_object.close()


def write_sam_line(sam_line=None, output_object=None):
    """Convert sam_dict to tab separated text and write out"""
    output_values = list(sam_line.values())
    sam1 = '\t'.join(output_values[:-1])
    sam2 = '\t'.join(output_values[-1])
    formatted_read = f'{sam1}\t{sam2}\n'
    output_object.write(formatted_read)


def write_bam_line(sam_line=None, output_object=None):
    """Convert sam line to bam line and write out"""
    output_values = list(sam_line.values())
    sam1 = '\t'.join(output_values[:-1])
    sam2 = '\t'.join(output_values[-1])
    formatted_read = f'{sam1}\t{sam2}'
    bam_line = pysam.AlignedSegment.fromstring(formatted_read, output_object.header)
    output_object.write(bam_line)
