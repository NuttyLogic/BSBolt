#! /usr/env python3

import gzip


class OpenAln:
    """ Simple class to simplify iterating through .aln files. The script yields a tuple for every thre lines
    in a fastq file
    ------------------------------------------------------------------------------------
    input: path to fastq
    returns: fastq iteration object"""

    def __init__(self, fastq=None):
        if fastq.endswith(".gz"):
            self.f = gzip.open(fastq, 'rb')
        else:
            self.f = open(fastq, 'r')

    def __iter__(self):
        with self.f as fastq:
            aln_start = False
            while True:
                line1 = fastq.readline()
                if not line1:
                    break
                if not aln_start:
                    if '##Header End' in line1:
                        aln_start = True
                    continue
                line1 = self.process_line(line1)
                line2 = self.process_line(fastq.readline())
                line3 = self.process_line(fastq.readline())
                yield (line1, line2, line3)

    @staticmethod
    def process_line(line):
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '')
        else:
            return line.replace('\n', '')
