import argparse
from BSB_CallMethylation.ProcessContigs import ProcessContigs
from BSB_Utils.BSB_UtilityFunctions import check_python_version


parser = argparse.ArgumentParser(description='BSBolt Module to Call Methylation Values from Sorted BAM Files')

parser.add_argument('-I', type=str, required=True,
                    help='Input BAM, input file must be in BAM format')
parser.add_argument('-DB', type=str, required=True, help='Path to index directory')
parser.add_argument('-O', type=str, required=True, help='Output prefix')
parser.add_argument('-remove-ccgg', action="store_true", default=False, help='Remove methylation calls in ccgg sites,'
                                                                             'default=False')
parser.add_argument('-verbose', action="store_true", default=False, help='Verbose Output, default=False')
parser.add_argument('-text', action="store_true", default=False, help='Output plain text files, default=False')
parser.add_argument('-remove-sx', action="store_false", default=True, help='Remove methylation calls from reads marked '
                                                                           'as incompletely by BSSeeker-Align, default='
                                                                           'True')
parser.add_argument('-ignore-overlap', action="store_true", default=False, help='Only consider higher quality base '
                                                                                'when paired end reads overlap, '
                                                                                'default=False')
parser.add_argument('-max', type=int, default=8000, help='Max read depth to call methylation')
parser.add_argument('-min', type=int, default=10, help='Minimum read depth required to report methylation site')
parser.add_argument('-t', type=int, default=1, help='Number of threads to use when calling methylation values')


if __name__ == "__main__":
    check_python_version()
    arguments = parser.parse_args()
    processed_contigs = ProcessContigs(input_file=arguments.I,
                                       genome_database=arguments.DB,
                                       output_prefix=arguments.O,
                                       remove_ccgg=arguments.remove_ccgg,
                                       remove_sx_reads=arguments.remove_sx,
                                       text_output=arguments.text,
                                       min_read_depth=arguments.min,
                                       threads=arguments.t,
                                       verbose=arguments.verbose)
