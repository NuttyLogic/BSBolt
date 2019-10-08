#! /usr/env python3

import gzip
import io


class OpenCGmap:

    def __init__(self, cgmap=None):
        if cgmap.endswith('.gz'):
            self.f = io.BufferedReader(gzip.open(cgmap, 'rb'))
        else:
            self.f = open(cgmap, 'r')

    def __iter__(self):
        with self.f as cg:
            while True:
                line = cg.readline()
                if not line:
                    break
                yield self.process_line(line)

    @staticmethod
    def process_line(line):
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '').split('\t')
        else:
            return line.replace('\n', '').split('\t')
