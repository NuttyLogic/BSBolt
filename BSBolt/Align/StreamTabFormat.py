#!/bin/python


import subprocess
import os


class StreamTab:
    """Launch tab format conversion and output
        Keyword Arguments:
           fastq1: fastq1
           fastq2: fastq2
    """

    def __init__(self, fastq1=None, fastq2=None, unstranded=False, no_conversion=False):
        assert isinstance(fastq1, str)
        self.fastq1 = fastq1
        if fastq2:
            assert isinstance(fastq2, str)
        self.fastq2 = fastq2
        self.python_path = self.get_python_version
        self.unstranded = unstranded
        self.no_conversion = no_conversion
        self.tab_format_path = self.get_tab_format_path
        self.tab_conversion = self.stream_tab

    @property
    def stream_tab(self):
        conversion_command: list = self.format_conversion_command
        return subprocess.Popen(conversion_command, stdout=subprocess.PIPE, universal_newlines=True)

    @property
    def get_python_version(self):
        python_path = 'python3'
        try:
            subprocess.Popen([python_path, '--version'], stdout=subprocess.PIPE)
        except FileNotFoundError:
            python_path = 'python'
        return python_path

    @property
    def get_tab_format_path(self):
        current_path = os.path.realpath(__file__).split('/')
        tab_format_path_list = current_path[0:-2]
        tab_format_path_list.extend(['Utils', 'TabFormatConversion.py'])
        return '/'.join(tab_format_path_list)

    @property
    def format_conversion_command(self):
        conversion_command = [self.python_path, self.tab_format_path, '-fq1', self.fastq1]
        if self.fastq2:
            conversion_command.extend(['-fq2', self.fastq2])
        if self.unstranded:
            conversion_command.append('-u')
        if self.no_conversion:
            conversion_command.append('-nc')
        return conversion_command
