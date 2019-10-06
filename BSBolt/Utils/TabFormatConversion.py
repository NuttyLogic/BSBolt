#! /usr/bin/env python3

import argparse
from BSBolt.Utils.FastqIterator import OpenFastq
from BSBolt.Utils.UtilityFunctions import check_python_version

check_python_version()


class ConvertFastq:
    """
    Class to launch external bash fastq conversion streaming. Simply converts input fastq to tab5 or tab6 format
    (tab5, 1 input files; tab6, 2 input files).
    Keyword Arguments:
        fastq1 (str): path to first input fastq file
        fastq2 (str): path to second input fastq file
        unstranded (bool): convert reverse complement of watson / crick strands
    Attributes:
        self.fastq1 (str): path to first input fastq file
        self.fastq2 (str): path to second input fastq file
        self.unstranded (bool): convert complement of watson / crick strands
    """

    def __init__(self, fastq1=None, fastq2=None, unstranded=False, no_conversion=False):
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.unstranded = unstranded
        self.conversion_rules = (('C', 'T'), ('G', 'A'))
        if no_conversion:
            self.conversion_rules = (('C', 'C'), ('G', 'G'))

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
            tab_line = self.tab_conversion(line, self.conversion_rules[0], 1)
            print(tab_line)
            if self.unstranded:
                tab_line2 = self.tab_conversion(line, self.conversion_rules[1], 2)
                print(tab_line2)

    def pipe_tab6(self):
        # tab five format is read_name\tsequence\tquality\tread_name2\tsequence2\tquality2\n
        fastq_iterator1 = OpenFastq(self.fastq1)
        fastq_iterator2 = OpenFastq(self.fastq2)
        for line1, line2 in zip(fastq_iterator1, fastq_iterator2):
            tab_line1 = self.tab_conversion(line1, self.conversion_rules[0], 1)
            tab_line2 = self.tab_conversion(line2, self.conversion_rules[1], 1)
            print(f'{tab_line1}\t{tab_line2}')
            if self.unstranded:
                tab_line3 = self.tab_conversion(line1, self.conversion_rules[1], 2)
                tab_line4 = self.tab_conversion(line2, self.conversion_rules[0], 2)
                print(f'{tab_line3}\t{tab_line4}')

    @staticmethod
    def tab_conversion(fastq_line, replacement_rule, read_conversion):
        # strip fastq specific formatting
        sequence = fastq_line[1].upper()
        name = fastq_line[0].replace('@', '').split('/')[0]
        name = f'{name.split(" ")[0]}_BSBolt_{read_conversion}_{replacement_rule[0]}2{replacement_rule[1]}_{sequence}'
        quality = fastq_line[3]
        if replacement_rule:
            sequence = sequence.replace(replacement_rule[0], replacement_rule[1])
        return f'{name}\t{sequence}\t{quality}'


# simple argparse to launch externally
parser = argparse.ArgumentParser(description='Iterate through fastq files, while yielding in-silico '
                                             'bisulfite converted reads')

parser.add_argument('-fq1', type=str, default=None)
parser.add_argument('-fq2', type=str, default=None)
parser.add_argument('-u', action='store_true', default=False)
parser.add_argument('-nc', action='store_true', default=False)

arguments = parser.parse_args()

if __name__ == "__main__":
    fastq_conversion = ConvertFastq(fastq1=arguments.fq1,
                                    fastq2=arguments.fq2,
                                    unstranded=arguments.u,
                                    no_conversion=arguments.nc)
    fastq_conversion.convert_reads()
