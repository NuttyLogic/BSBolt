#! /usr/env python3

import gzip
import io
from typing import Tuple, Union


class OpenAln:
    """ Simple class to simplify iterating through .aln files.
    ------------------------------------------------------------------------------------
    input: path to fastq
    returns: fastq iteration object"""

    def __init__(self, aln: str = None):
        if aln.endswith(".gz"):
            self.f = io.BufferedReader(gzip.open(aln, 'rb'))
        else:
            self.f = open(aln, 'r')

    def __iter__(self) -> Tuple[str, str, str]:
        with self.f as aln:
            aln_start = False
            while True:
                line1 = aln.readline()
                if not line1:
                    break
                if not aln_start:
                    if '##Header End' in line1:
                        aln_start = True
                    continue
                line1 = self.process_line(line1)
                line2 = self.process_line(aln.readline())
                line3 = self.process_line(aln.readline())
                yield line1, line2, line3

    @staticmethod
    def process_line(line: Union[str, bytes]) -> str:
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '')
        else:
            return line.replace('\n', '')
