#! /usr/env python3

from BSB.BSB_Align.StreamTabFormat import StreamTab
from BSB.BSB_Utils.SamIterator import OpenSam


class TabSamIterator:
    """Wrapper to iterate through fastq and sam lines in unison. """

    def __init__(self, fastq1=None, fastq2=None, sam_tuple=None):
        assert isinstance(sam_tuple, tuple)
        self.paired_end = False
        if fastq2:
            self.paired_end = True
        self.tab_iterator = StreamTab(fastq1=fastq1, fastq2=fastq2)
        self.sam_iterators = self.get_sam_iterator(sam_tuple=sam_tuple)
        self.iteration_order = ('W_C2T', 'C_C2T', 'W_G2A', 'C_G2A')

    def get_sam_iterator(self, sam_tuple):
        sam_iterators = []
        for sam_file in sam_tuple:
            sam_iterators.append(OpenSam(sam_file=sam_file, paired_end=self.paired_end))
        return sam_iterators

    def __iter__(self):
        for line in zip(self.tab_iterator, *self.sam_iterators):
            read_count = 1
            for read_id, read_sequence in line[0].items():
                if self.paired_end:
                    read_id = read_id.split('/')[0]
                    read_id = read_id.split(' ')[0]
                    read_id = f'{read_id}_{read_count}'
                output_dict = {'read_sequence': read_sequence}
                for reference_type, read_dictionary in zip(self.iteration_order, line[1:]):
                    output_dict.update({reference_type: read_dictionary[read_id]})
                yield output_dict
                read_count += 1
