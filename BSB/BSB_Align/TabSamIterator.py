#! /usr/env python3

from BSB.BSB_Utils.SamIterator import OpenSam


class TabSamIterator:
    """Wrapper to iterate through fastq and sam lines in unison. """

    def __init__(self, fastq1=None, fastq2=None, sam_tuple=None):
        assert isinstance(sam_tuple, tuple)
        self.paired_end = False
        if fastq2:
            self.paired_end = True
        self.sam_iterators = self.get_sam_iterator(sam_tuple=sam_tuple)
        self.iteration_order = ('W_C2T', 'C_C2T', 'W_G2A', 'C_G2A')

    def get_sam_iterator(self, sam_tuple):
        sam_iterators = []
        for sam_file in sam_tuple:
            sam_iterators.append(OpenSam(sam_file=sam_file, paired_end=self.paired_end))
        return sam_iterators

    @staticmethod
    def process_read_id(read_identifier):
        read_id = read_identifier.split('_BSBolt_')[0]
        read_id = read_id.split('/')[0]
        read_id = read_id.split(' ')[0]
        return read_id

    @staticmethod
    def get_ordered_reads(sam_reads, max_read_number):
        ordered_reads = []
        for read_index in range(max_read_number):
            ordered_sam_reads = []
            for sam_index, sam_grouping in enumerate(sam_reads):
                try:
                    ordered_sam_reads.append(sam_grouping[read_index])
                except IndexError:
                    ordered_sam_reads.append(sam_grouping[0])
            ordered_reads.append(ordered_sam_reads)
        return ordered_reads

    def __iter__(self):
        for line in zip(*self.sam_iterators):
            max_read_number = max([len(reference_read) for reference_read in line])
            ordered_reads = self.get_ordered_reads(line, max_read_number)
            original_sequences = [None, None]
            for sam_reads in ordered_reads:
                output_dicts = [dict(read_sequence=None), dict(read_sequence=None)]
                for reference_type, read_dictionary in zip(self.iteration_order, sam_reads):
                    for index, sam_info in enumerate(read_dictionary.values()):
                        if 'original_sequence' in sam_info:
                            original_sequences[index] = str(sam_info['original_sequence'])
                            del sam_info['original_sequence']
                        if not output_dicts[index]['read_sequence']:
                            output_dicts[index]['read_sequence'] = original_sequences[index]
                        output_dicts[index].update({reference_type: sam_info})
                yield output_dicts[0]
                if len(output_dicts[1]) > 1:
                    yield output_dicts[1]
