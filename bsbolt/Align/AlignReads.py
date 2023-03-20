import os
import subprocess
import sys
from typing import List
from bsbolt.Utils.UtilityFunctions import get_external_paths


class BisulfiteAlignmentError(Exception):
    """Error in alignment"""
    pass


class AlignmentCompressionError(Exception):
    """Error in read compression"""
    pass


class BisulfiteAlignmentAndProcessing:
    """ Read alignment as processing. Handles reads from BWA-MEM2 and outputs to BAM file. Also aggregates alignment
    stats.

    Params:

    * *alignment_commands (list)*: bwa alignment commands
    * *output (str)*: output prefix
    * *output_threads (int)*: number of threads available for bam output
    * *output_to_stdout (boold)*: output alignments to stdout

    Attributes:

    * self.mapping_statistics (dict)*: alignment run statistics
    """

    def __init__(self, alignment_commands: List[str], output: str = None, output_threads: int = 1,
                 output_to_stdout: bool = False):
        self.alignment_commands = alignment_commands
        self.output = output
        self.output_threads = output_threads
        self.output_to_stdout = output_to_stdout
        self.mapping_statistics = dict(TotalReads=0, TotalAlignments=0, BSAmbiguous=0, C_C2T=0, C_G2A=0,
                                       W_C2T=0, W_G2A=0, Unaligned=0)

    def align_reads(self):
        """ Launch bwa alignment. Pipe output to BAM file
        """
        _, _, stream_bam = get_external_paths()
        if '/' in self.output:
            assert os.path.exists('/'.join(self.output.split('/')[0:-1])), f"output path {self.output} not valid"
        bam_compression = None
        if not self.output_to_stdout:
            alignment_run = subprocess.Popen(self.alignment_commands,
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE,
                                             universal_newlines=True)
            bam_compression = subprocess.Popen([stream_bam, '-@', str(self.output_threads), '-o', f'{self.output}.bam'],
                                                stdin=alignment_run.stdout)
        else:
            alignment_run = subprocess.Popen(self.alignment_commands,
                                             stderr=subprocess.PIPE,
                                             universal_newlines=True)
        # watch alignment progress, output stderr and collect alignment stats
        while True:
            # Show intermediate steps of alignment
            alignment_info = alignment_run.stderr.readline()
            if alignment_info:
                if alignment_info[0:7] == 'BSStat ':
                    category, count = alignment_info.replace('BSStat ', '').split(': ')
                    self.mapping_statistics[category] += int(count)
                    print(alignment_info.replace('BSStat ', '').strip(), file=sys.stderr)
                else:
                    print(alignment_info.strip(), file=sys.stderr)
            if bam_compression is not None:
                if bam_compression.returncode:
                    print(bam_compression.returncode, file=sys.stderr)
                    raise AlignmentCompressionError
            if alignment_run.returncode:
                print(alignment_run.returncode, file=sys.stderr)
                raise BisulfiteAlignmentError
            if bam_compression is not None:
                if alignment_run.poll() is not None and bam_compression.poll() is not None:
                    alignment_run.stdout.close()
                    alignment_run.stderr.close()
                    break
            else:
                if alignment_run.poll() is not None:
                    alignment_run.stderr.close()
                    break
