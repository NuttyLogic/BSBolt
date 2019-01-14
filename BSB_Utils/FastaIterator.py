#! /usr/env python3

import gzip


class OpenFasta:
    """ Simple class to simplify iterating through fasta
    ------------------------------------------------------------------------------------
    input: path to fastq
    returns: fastq iteration object"""

    def __init__(self, fasta=None):
        if fasta.endswith(".gz"):
            self.f = gzip.open(fasta, 'rb')
        else:
            self.f = open(fasta, 'r')

    def __iter__(self):
        with self.f as fasta:
            while True:
                line = fasta.readline()
                if not line:
                    break
                if '>' in line:
                    yield True, self.process_line(line)
                else:
                    yield False, self.process_line(line)

    @staticmethod
    def process_line(line):
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '')
        else:
            return line.replace('\n', '')
