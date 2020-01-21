import argparse

from BSBolt.Utils.UtilityFunctions import get_external_paths

bt2_path, art_path = get_external_paths()

parser = argparse.ArgumentParser(description='BiSulfite Bolt, A bisulfite sequencing processing platform.',
                                 usage='BSBolt Module {Module Arguments}')

subparsers = parser.add_subparsers(description='BSBolt Modules, Please Invoke BSBolt Modules for Help\n '
                                               'Documentation at bsbolt.readthedocs.io',
                                   metavar='Index, Align, CallMethylation, AggregateMatrix, Simulate, Impute',
                                   dest='subparser_name')

align_parser = subparsers.add_parser('Align', help='Alignment Module')
index_parser = subparsers.add_parser('Index', help='Index Generation Module')
call_meth_parser = subparsers.add_parser('CallMethylation', help='Methylation Calling Module')
matrix_parser = subparsers.add_parser('AggregateMatrix', help='CGmap Matrix Aggregation Module')
sim_parser = subparsers.add_parser('Simulate', help='BSBolt Illumina Read Simulation Module')
imputation_parser = subparsers.add_parser('Impute', help='kNN Imputation Module')

# Add Alignment Parser Commands

align_parser.add_argument('-F1', type=str, default=None, help='Path to fastq 1', required=True)
align_parser.add_argument('-F2', type=str, default=None, help='Path to fastq 2')
align_parser.add_argument('-NC', action='store_true', default=False, help='Aligned unconverted bisulfite reads')
align_parser.add_argument('-D', action="store_true", default=False, help='Library directional, in silico modification '
                                                                         'of sequencing bases will only account for '
                                                                         'bisulfite converted DNA and not the PCR '
                                                                         'products of converted DNA, default=False')
align_parser.add_argument('-O', type=str, default=None, help='Path to Output Prefix', required=True)
align_parser.add_argument('-DB', type=str, default=None, help='Path to BSBolt Database', required=True)
align_parser.add_argument('-M', type=int, default=4, help='Read mismatch threshold, reads where the total number of '
                                                          'observed mismatches greater than the threshold will be '
                                                          'discarded')
align_parser.add_argument('-S', action="store_true", default=False, help='Position sort output .bam, default=False')
align_parser.add_argument('-BT2-local', action="store_true", default=False,
                          help='Bowtie2; local alignment, default end-to-end')
align_parser.add_argument('-BT2-D', type=int, default=40, help='Bowtie2; number of consecutive seed extension attempts '
                                                               'that can fail before Bowtie2 move on')
align_parser.add_argument('-BT2-k', type=int, default=5, help='Bowtie2; alignment search limit')
align_parser.add_argument('-BT2-p', type=int, default=2, help='Bowtie2; Number of threads for Bowtie2 to use')
align_parser.add_argument('-BT2-L', type=int, default=20, help='Bowtie2; Length of subseeds during alignment')
align_parser.add_argument('-BT2-score-min', type=str, default='L,-0.6,-0.6', help='Bowtie2; scoring function')
align_parser.add_argument('-BT2-I', type=int, default=0, help='Bowtie2; minimum fragment length '
                                                              'for a valid paired-end alignment')
align_parser.add_argument('-BT2-X', type=int, default=500, help='Bowtie2; maximum fragment length '
                                                                'for a valid paired-end alignment')
align_parser.add_argument('-discordant', action='store_true', help='Report discordant and mixed reads, default=False')

# Add Index Parser Commands

index_parser.add_argument('-G', type=str, required=True,
                          help='Path to reference genome fasta file, fasta file should contain all contigs')
index_parser.add_argument('-DB', type=str, required=True,
                          help='Path to index directory, will create directory if folder does not exist')
index_parser.add_argument('-MR', type=str, default=None, help='Path to bed file of mappable regions. Index will be '
                                                              'built using using masked contig sequence')
index_parser.add_argument('-BT2-p', type=int, default=2, help='Bowtie2; Number of threads')
index_parser.add_argument('-rrbs', action="store_true", default=False, help='Generate a Reduced Representative'
                                                                            ' Bisulfite Sequencing (RRBS) index')
index_parser.add_argument('-rrbs-cut-format', default='C-CGG',
                          help='Cut format to use for generation of RRBS database, '
                               'default= C-CGG (MSPI), input multiple enzymes as a '
                               'comma separate string, C-CGG,C-CGG,...')
index_parser.add_argument('-rrbs-lower', type=int, default=40, help='Lower bound fragment size to consider RRBS index'
                                                                    'generation, default = 40')
index_parser.add_argument('-rrbs-upper', type=int, default=500, help='Upper bound fragment size to consider RRBS index'
                                                                     'generation, default = 500')

# Add Methylation Calling Parser Arguments

call_meth_parser.add_argument('-I', type=str, required=True,
                              help='Input BAM, input file must be in BAM format with index file')
call_meth_parser.add_argument('-DB', type=str, required=True, help='Path to index directory')
call_meth_parser.add_argument('-O', type=str, required=True, help='Output prefix')
call_meth_parser.add_argument('-remove-ccgg', action="store_true", default=False,
                              help='Remove methylation calls in ccgg sites,'
                                   'default=False')
call_meth_parser.add_argument('-verbose', action="store_true", default=False, help='Verbose Output, default=False')
call_meth_parser.add_argument('-text', action="store_true", default=False,
                              help='Output plain text files, default=False')
call_meth_parser.add_argument('-ignore-overlap', action="store_true", default=False,
                              help='Only consider higher quality base '
                                   'when paired end reads overlap, '
                                   'default=False')
call_meth_parser.add_argument('-max', type=int, default=8000, help='Max read depth to call methylation')
call_meth_parser.add_argument('-min', type=int, default=10,
                              help='Minimum read depth required to report methylation site')
call_meth_parser.add_argument('-t', type=int, default=1,
                              help='Number of threads to use when calling methylation values')
call_meth_parser.add_argument('-min-qual', type=int, default=0, help='Minimum base quality for a base to considered for'
                                                                     'methylation calling, default=0')
call_meth_parser.add_argument('-CG', action="store_true", default=False, help='Only output CpG sites in CGmap file')
call_meth_parser.add_argument('-ATCG', action="store_true", default=False, help='Output ATCGmap file')
call_meth_parser.add_argument('-IO', action="store_true", default=False, help='Ignore orphans during methylation call')

# Add Matrix Aggregation Parser Args

matrix_parser.add_argument('-F', type=lambda file: [file_path for file_path in file.split(',')], required=True,
                           help='Comma separated list of CGmap file paths, or '
                                'path to text file with list of line separated '
                                'CGmap file paths')
matrix_parser.add_argument('-S', type=lambda sample_labels: [sample for sample in sample_labels.split(',')],
                           default=None,
                           help='Comma separated list of samples labels. '
                                'If sample labels are not provided sample labels '
                                'are extracted from CGmap file paths. '
                                'Can also pass path to txt for line separated sample '
                                'labels.')
matrix_parser.add_argument('-min-coverage', type=int, default=10, help='Minimum site read depth coverage for a '
                                                                       'site to be included in the aggregate matrix')
matrix_parser.add_argument('-min-sample', type=float, default=0.80,
                           help='Proportion of samples that must have a valid site '
                                '(above minimum coverage threshold), for a site to be'
                                'included in the aggregate matrix.')
matrix_parser.add_argument('-O', type=str, default=None, required=True, help='Aggregate matrix output path')
matrix_parser.add_argument('-CG', action="store_true", default=False, help='Only output CG sites')
matrix_parser.add_argument('-verbose', action="store_true", default=False, help='Verbose aggregation')
matrix_parser.add_argument('-t', type=int, default=1, help='Number of threads to use when assembling matrix')
matrix_parser.add_argument('-count', action='store_true', default=False,
                           help='Output a count matrix with count of methylated cytosines and '
                                'total observed cytosines, Sample_meth_cytosines   Sample_total_cytosines')


# Add Simulation Parser Args

sim_parser.add_argument('-G', type=str, required=True,
                        help='Path for reference genome fasta file, fasta file should contain all contigs')
sim_parser.add_argument('-O', type=str, required=True, help='Output prefix')
sim_parser.add_argument('-PE', default=False, action='store_true', help='Simulate Paired End Reads, default Single End')
sim_parser.add_argument('-RL', type=int, default=125, help='Simulated Read Length')
sim_parser.add_argument('-RD', type=int, default=20, help='Simulated Read Depth')
sim_parser.add_argument('-U', default=False, action='store_true',
                        help='Simulate Undirectional Reads, default=Directional')
sim_parser.add_argument('-RC', default=None, help='Path to CGmap file to generate simulation reference profile')
sim_parser.add_argument('-RO', default=None, help='Methylation reference output directory, default = output path')
sim_parser.add_argument('-BR', default=None, help='Path to previously generate BSBolt simulation reference')
sim_parser.add_argument('-IR1', type=float, default=0.001, help='Read 1 insertion rate')
sim_parser.add_argument('-IR2', type=float, default=0.001, help='Read 2 insertion rate')
sim_parser.add_argument('-DR1', type=float, default=0.001, help='Read 1 deletion rate')
sim_parser.add_argument('-DR2', type=float, default=0.001, help='Read 2 deletion rate')
sim_parser.add_argument('-NF', type=int, default=0, help='Cutoff threshold for a read with gaps, -, or ambiguous bases, '
                                                         'N. Reads below the threshold will not be output')
sim_parser.add_argument('-M', type=int, default=400, help='Mean paired end fragment size')
sim_parser.add_argument('-SM', type=int, default=50, help='Paired end fragment length distribution standard deviation')
sim_parser.add_argument('-SS', type=str, default='HS25', help="Sequencing system, default HS25 (HiSeq 2500)\n "
                                                              "HS10 - HiSeq 1000 (100bp) \n\
                                                               HS20 - HiSeq 2000 (100bp)\n\
                                                               HS25 - HiSeq 2500 (125bp, 150bp)\n\
                                                               HSXn - HiSeqX PCR free (150bp)\n\
                                                               HSXt - HiSeqX TruSeq (150bp)\n\
                                                               MinS - MiniSeq TruSeq (50bp)\n\
                                                               MSv1 - MiSeq v1 (250bp)\n\
                                                               MSv3 - MiSeq v3 (250bp)\n\
                                                               NS50 - NextSeq500 v2 (75bp)\n")
sim_parser.add_argument('-Q1', type=str, default=None, help='Optional read 1 quality profile, generated using ART '
                                                            'Illumina read profiler')
sim_parser.add_argument('-Q2', type=str, default=None, help='Path to read 2 quality profile')

# Add Imputation Parser Args

imputation_parser.add_argument('-M', type=str, required=True, help='Path to BSBolt matrix file')
imputation_parser.add_argument('-B', type=int, default=0, help='Imputation sample batch size kNN imputation, by default'
                                                               ' the all of the samples will be processed as a single '
                                                               'batch')
imputation_parser.add_argument('-W', type=int, default=3000000, help='Sliding window size for imputation')
imputation_parser.add_argument('-k', type=int, default=5, help='Number of neighbors to use for imputation, default = 5')
imputation_parser.add_argument('-t', type=int, default=1, help='Number of threads available for imputation')
imputation_parser.add_argument('-verbose', action='store_true', default=False, help='Verbose output')
imputation_parser.add_argument('-O', type=str, default='', help='Output path for imputed matrix')
imputation_parser.add_argument('-R', action='store_true', default=False, help='Randomize batches')
