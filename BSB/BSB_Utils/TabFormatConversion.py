#! /usr/bin/env python3

import argparse
from FastqIterator import OpenFastq
from UtilityFunctions import check_python_version

check_python_version()


class ConvertFastq:
    """
    Class to launch external bash fastq conversion streaming. Simply converts input fastq to tab5 or tab6 format
    (tab5, 1 input files; tab6, 2 input files).
    Keyword Arguments:
        fastq1 (str): path to first input fastq file
        fastq2 (str): path to second input fastq file
        replacement_base1 (str): base to substitute
        replacement_base2 (str): base to substitute
    Attributes:
        self.fastq1 (str): path to first input fastq file
        self.fastq2 (str): path to second input fastq file
        self.replacement_base1 (str): base to substitute
        self.replacement_base2 (str): base to substitute
    """

    def __init__(self, fastq1=None, fastq2=None, replacement_base1=None, replacement_base2=None):
        assert isinstance(fastq1, str)
        self.fastq1 = fastq1
        if fastq2:
            assert isinstance(fastq2, str)
        self.fastq2 = fastq2
        if replacement_base1:
            assert isinstance(replacement_base1, str)
            replacement_base1 = self.get_replacemet_rule(str(replacement_base1))
        self.replacement_base1 = replacement_base1
        if replacement_base2:
            assert isinstance(replacement_base2, str)
            replacement_base2 = self.get_replacemet_rule(str(replacement_base2))
        self.replacement_base2 = replacement_base2

    def convert_reads(self):
        # output tab5 if self.fastq2 not true
        if not self.fastq2:
            self.pipe_tab5()
        else:
            self.pipe_tab6()

    def pipe_tab5(self):
        # tab five format is read_name\tsequence\tquality\n
        fastq_iterator1 = OpenFastq(self.fastq1)
        for line in fastq_iterator1:
            tab_line = self.tab_conversion(line, self.replacement_base1)
            print(f'{tab_line}')

    def pipe_tab6(self):
        # tab five format is read_name\tsequence\tquality\tread_name2\tsequence2\tquality2\n
        fastq_iterator1 = OpenFastq(self.fastq1)
        fastq_iterator2 = OpenFastq(self.fastq2)
        for line1, line2 in zip(fastq_iterator1, fastq_iterator2):
            tab_line1 = self.tab_conversion(line1, self.replacement_base1)
            tab_line2 = self.tab_conversion(line2, self.replacement_base2)
            print(f'{tab_line1}\t{tab_line2}')

    @staticmethod
    def tab_conversion(fastq_line, replacement_rule):
        # strip fastq specific formatting
        sequence = fastq_line[1]
        name = fastq_line[0].replace('@', '').split('/')[0]
        name = f'{name.split(" ")[0]}_BSBolt_{sequence}'
        quality = fastq_line[3]
        if replacement_rule:
            sequence = sequence.replace(replacement_rule[0], replacement_rule[1].lower())
        return f'{name}\t{sequence}\t{quality}'

    @staticmethod
    def get_replacemet_rule(replacement_base):
        # only replace Cytosine or Guanine for purposes of methylation alignment
        if replacement_base != 'C' and replacement_base != 'G':
            raise TypeError('Replacement Base Must by Cytosine or Guanine')
        replacement_rule = ('C', 'T')
        if replacement_base == 'G':
            replacement_rule = ('G', 'A')
        return replacement_rule


# simple argparse to launch externally
parser = argparse.ArgumentParser(description='Iterate through fastq files, while yielding in-silico '
                                             'bisulfite converted reads')

parser.add_argument('-fq1', type=str, default=None)
parser.add_argument('-fq2', type=str, default=None)
parser.add_argument('-r1', type=str, default=None)
parser.add_argument('-r2', type=str, default=None)

arguments = parser.parse_args()


if __name__ == "__main__":
    fastq_conversion = ConvertFastq(fastq1=arguments.fq1,
                                    fastq2=arguments.fq2,
                                    replacement_base1=arguments.r1,
                                    replacement_base2=arguments.r2)
    fastq_conversion.convert_reads()
