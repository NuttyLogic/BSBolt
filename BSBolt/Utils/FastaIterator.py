#! /usr/env python3

import gzip
import io
from typing import Tuple, Union


class OpenFasta:
    """ Simple class to simplify iterating through fasta
    ------------------------------------------------------------------------------------
    KeywordArguments:
        fasta: str = path to fastq
    """

    def __init__(self, fasta: str = None):
        if fasta.endswith(".gz"):
            self.f = io.BufferedReader(gzip.open(fasta, 'rb'))
        else:
            self.f = open(fasta, 'r')

    def __iter__(self) -> Tuple[bool, str]:
        with self.f as fasta:
            while True:
                line = fasta.readline()
                if not line:
                    break
                processed_line = self.process_line(line)
                if '>' in processed_line:
                    yield True, processed_line
                else:
                    yield False, processed_line

    @staticmethod
    def process_line(line: Union[str, bytes]) -> str:
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '')
        else:
            return line.replace('\n', '')
