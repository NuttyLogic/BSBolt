#! /usr/env python3

import gzip
import io


class OpenFastq:
    """ Simple class to simplify iterating through fastq files. The script yields a tuple for every four lines
    in a fastq file
    ------------------------------------------------------------------------------------
    input: path to fastq
    returns: fastq iteration object"""

    def __init__(self, fastq=None):
        if fastq.endswith(".gz"):
            self.f = io.BufferedReader(gzip.open(fastq, 'rb'))
        else:
            self.f = open(fastq, 'r')

    def __iter__(self):
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
    def process_line(line):
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '')
        else:
            return line.replace('\n', '')
