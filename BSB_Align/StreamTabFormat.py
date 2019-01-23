import subprocess
import os
from distutils.version import StrictVersion


class StreamTab:
    """Launch tab format conversion and output
        Keyword Arguments:
           fastq1: fastq1
           fastq2: fastq2
           replacement_base1: base to replace
           replacement_base2: base to replace 2
    """

    def __init__(self, fastq1=None, fastq2=None, replacement_base1=None, replacement_base2=None):
        assert isinstance(fastq1, str)
        self.fastq1 = fastq1
        if fastq2:
            assert isinstance(fastq2, str)
        self.fastq2 = fastq2
        if replacement_base1:
            assert isinstance(replacement_base1, str)
        self.replacement_base1 = replacement_base1
        if replacement_base2:
            assert isinstance(replacement_base2, str)
        self.replacement_base2 = replacement_base2
        self.python_path = self.get_python_version
        self.tab_format_path = self.get_tab_format_path
        self.tab_conversion = self.stream_tab

    @property
    def stream_tab(self):
        conversion_command: list = self.format_conversion_command
        return subprocess.Popen(conversion_command, stdout=subprocess.PIPE, universal_newlines=True)

    def __iter__(self):
        for tab_line in iter(self.tab_conversion.stdout.readline, ''):
            if tab_line != '\n':
                yield self.convert_tab_format(tab_line)
        self.tab_conversion.stdout.close()
        return_code = self.tab_conversion.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, 'Pipe error')

    @property
    def get_python_version(self):
        python_path = 'python3'
        try:
            version_command = subprocess.Popen([python_path, '--version'], stdout=subprocess.PIPE)
        except FileNotFoundError:
            python_path = 'python'
        return python_path

    @property
    def get_tab_format_path(self):
        current_path = os.path.realpath(__file__).split('/')
        tab_format_path_list = current_path[0:-2]
        tab_format_path_list.extend(['BSB_Utils', 'TabFormatConversion.py'])
        return '/'.join(tab_format_path_list)

    @property
    def format_conversion_command(self):
        conversion_command = [self.python_path, self.tab_format_path, '-fq1', self.fastq1]
        if self.replacement_base1:
            conversion_command.extend(['-r1', self.replacement_base1])
        if self.replacement_base2:
            conversion_command.extend(['-r2', self.replacement_base2])
        if self.fastq2:
            conversion_command.extend(['-fq2', self.fastq2])
        return conversion_command

    def convert_tab_format(self, tab_line):
        tab_split: list = tab_line.replace('\n', '').split('\t')
        if self.fastq2:
            return {f'{tab_split[0]}/1': [tab_split[1], tab_split[2]], f'{tab_split[3]}/2': [tab_split[4], tab_split[5]]}
        return {tab_split[0]: [tab_split[1], tab_split[2]]}
