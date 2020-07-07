import os
import subprocess
from typing import List
from BSBolt.Utils.UtilityFunctions import get_external_paths


class BisulfiteAlignmentAndProcessing:
    """ Read alignment as processing. Handles reads from BWA-MEM2 and outputs to BAM file. Also aggregates alignment
    stats.

    Params:

    * *alignment_commands (list)*: bwa alignment commands
    * *output (str)*: output prefix

    Attributes:

    * self.mapping_statistics (dict)*: alignment run statistics
    """

    def __init__(self, alignment_commands: List[str], output=None):
        self.alignment_commands = alignment_commands
        self.output = output
        self.mapping_statistics = dict(TotalReads=0, TotalAlignments=0, BSAmbiguous=0, C_C2T=0, C_G2A=0,
                                       W_C2T=0, W_G2A=0, Unaligned=0)

    def align_reads(self):
        """ Launch bwa alignment. Pipe output to BAM file
        """
        _, _, stream_bam = get_external_paths()
        if '/' in self.output:
            assert os.path.exists('/'.join(self.output.split('/')[0:-1])), f"output path {self.output} not valid"
        alignment_run = subprocess.Popen(self.alignment_commands,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE,
                                         universal_newlines=True)
        bam_compression = subprocess.Popen([stream_bam, '-o', f'{self.output}.bam'],
                                           stdin=alignment_run.stdout)
        # watch alignment progress, output stderr and collect alignment stats
        while True:
            if bam_compression.returncode:
                break
            elif alignment_run.returncode:
                break
            elif alignment_run.poll() is not None and bam_compression.poll() is not None:
                alignment_run.stdout.close()
                alignment_run.stderr.close()
                break
            else:
                alignment_info = alignment_run.stderr.readline().strip()
                if alignment_info:
                    if alignment_info[0:7] == 'BSStat ':
                        category, count = alignment_info.replace('BSStat ', '').split(': ')
                        self.mapping_statistics[category] += int(count)
                    else:
                        print(alignment_info)
