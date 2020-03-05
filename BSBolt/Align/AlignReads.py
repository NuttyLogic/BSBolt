import subprocess
from typing import Any, Dict, List, Set
import pysam


class BisulfiteAlignmentAndProcessing:
    """
    """

    def __init__(self,  alignment_commands: List[str], output=None):
        self.alignment_commands = alignment_commands
        self.output = output
        self.mapping_statistics = dict(TotalReads=0, TotalAlignments=0, BSAmbiguous=0, C_C2T=0, C_G2A=0,
                                       W_C2T=0, W_G2A=0, Unaligned=0)

    def align_reads(self):
        """
        """
        alignment_run = subprocess.Popen(self.alignment_commands,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         universal_newlines=True)
        infile = pysam.AlignmentFile(alignment_run.stdout, 'r')
        sam_out = pysam.AlignmentFile(f'{self.output}.bam', 'wb', template=infile)
        for s in infile:
            sam_out.write(s)
        self.update_mapping_statistics(alignment_run.stderr)

    def update_mapping_statistics(self, log_pipe):
        """Mapping statistics are updated per read group based on the mapping status"""
        for line in iter(log_pipe.readline, ''):
            if line[0:7] == 'BSStat ':
                category, count = line.replace('BSStat ', '').split(': ')
                self.mapping_statistics[category] += int(count)
            else:
                if line:
                    print(line.strip())
