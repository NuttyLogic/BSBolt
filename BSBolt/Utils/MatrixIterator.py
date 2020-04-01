#! /usr/env python3

import gzip
import io
from typing import List, Tuple, Union
import numpy as np


class OpenMatrix:
    """ Simple class to simple class to iterate through BSBolt matrix file
    ------------------------------------------------------------------------------------
    input: path to matrix file
    returns: matrix iteration object"""

    def __init__(self, matrix: str = None):
        self.header = False
        if matrix.endswith(".gz"):
            self.f = io.BufferedReader(gzip.open(matrix, 'rb'))
        else:
            self.f = open(matrix, 'r')

    def __iter__(self):
        with self.f as matrix:
            while True:
                line = matrix.readline()
                if not line.strip():
                    break
                line = self.process_line(line)
                yield line

    def process_line(self, line) -> Tuple[str, Union[List[str], np.ndarray]]:
        converted_line = self.line_conversion(line)
        if not self.header:
            self.header = True
            return converted_line[0], converted_line[1:]
        return converted_line[0], np.asarray([self.convert_to_float(value) for value in converted_line[1:]])

    @staticmethod
    def line_conversion(line) -> List[str]:
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '').split('\t')
        else:
            return line.replace('\n', '').split('\t')

    @staticmethod
    def convert_to_float(value: str) -> Union[float, None]:
        try:
            return float(value)
        except ValueError:
            return np.nan
