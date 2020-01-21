#! /usr/env python3

import gzip
import io
from typing import List, Union


class OpenFastq:
    """ Simple class to simplify iterating through fastq files. The script yields a tuple for every four lines
    in a fastq file
    ------------------------------------------------------------------------------------
    KeywordArguments:
        fastq: str = path to fastq file
    """

    def __init__(self, fastq: str = None):
        if fastq.endswith(".gz"):
            self.f = io.BufferedReader(gzip.open(fastq, 'rb'))
        else:
            self.f = open(fastq, 'r')

    def __iter__(self) -> List[str]:
        with self.f as fastq:
            while True:
                line1 = fastq.readline()
                if not line1.strip():
                    break
                line1 = self.process_line(line1)
                line2 = self.process_line(fastq.readline())
                line3 = self.process_line(fastq.readline())
                line4 = self.process_line(fastq.readline())
                yield [line1, line2, line3, line4]

    @staticmethod
    def process_line(line: Union[bytes, str]):
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '')
        else:
            return line.replace('\n', '')
