#! /usr/env python3

import gzip
import io
import numpy as np


class OpenMatrix:
    """ Simple class to simple class to iterate through BSB matrix file
    ------------------------------------------------------------------------------------
    input: path to matrix file
    returns: matrix iteration object"""

    def __init__(self, matrix=None):
        self.header = None
        if matrix.endswith(".gz"):
            self.f = io.BufferedReader(gzip.open(matrix, 'rb'))
        else:
            self.f = open(matrix, 'r')

    def __iter__(self):
        with self.f as fastq:
            while True:
                line1 = fastq.readline()
                if not line1.strip():
                    break
                line1 = self.process_line(line1)
                yield line1

    def process_line(self, line):
        converted_line = self.line_conversion(line)
        if not self.header:
            self.header = 'Done'
            return converted_line[0], converted_line[1:]
        return converted_line[0], np.asarray(self.convert_to_float(value) for value in converted_line[1:])

    @staticmethod
    def line_conversion(line):
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '').split('\t')
        else:
            return line.replace('\n', '').split('\t')

    @staticmethod
    def convert_to_float(value):
        try:
            return float(value)
        except ValueError:
            return np.nan
