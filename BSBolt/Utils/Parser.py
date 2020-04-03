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
sort_parser = subparsers.add_parser('Sort', help='BAM Sort')
bam_index = subparsers.add_parser('BamIndex', help='BAM Index')
# Add Alignment Parser Commands

align_parser.add_argument('-F1', type=str, default=None, help='Path to fastq 1', required=True)
align_parser.add_argument('-F2', type=str, default=None, help='Path to fastq 2')
align_parser.add_argument('-UN', action="store_true", default=False,
                          help='Library Undirectional, Consider PCR products of bisulfite converted DNA',
                          required=False)
align_parser.add_argument('-O', type=str, default=None, help='Path to Output Prefix', required=True)
align_parser.add_argument('-DB', type=str, default=None, help='Path to BSBolt Database', required=True)
align_parser.add_argument('-t', type=int, default=1, help='Number of bwa threads [1]', required=False)
align_parser.add_argument('-k', type=int, default=19, help='Minimum seed length [19]', required=False)
align_parser.add_argument('-w', type=int, default=100, help='Band width for banded alignment [100]', required=False)
align_parser.add_argument('-d', type=int, default=100, help='off-diagonal X-dropoff [100]', required=False)
align_parser.add_argument('-r', type=float, default=1.5,
                          help='look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]', required=False)
align_parser.add_argument('-y', type=int, default=20,
                          help='seed occurrence for the 3rd round seeding [20]', required=False)
align_parser.add_argument('-c', type=int, default=500,
                          help='skip seeds with more than INT occurrences [500]', required=False)
align_parser.add_argument('-D', type=float, default=0.50,
                          help='drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]',
                          required=False)
align_parser.add_argument('-W', type=int, default=0, help='discard a chain if seeded bases shorter than INT [0]',
                          required=False)
align_parser.add_argument('-m', type=int, default=50,
                          help='perform at most INT rounds of mate rescues for each read [50]',
                          required=False)
align_parser.add_argument('-S', action='store_true', default=False, help='skip mate rescue', required=False)
align_parser.add_argument('-P', action='store_true', default=False,
                          help='skip pairing; mate rescue performed unless -S also in use', required=False)
align_parser.add_argument('-A', type=int, default=1,
                          help='score for a sequence match, which scales options -TdBOELU unless overridden [1]',
                          required=False)
align_parser.add_argument('-B', type=int, default=4,
                          help='penalty for a mismatch [4]',
                          required=False)
align_parser.add_argument('-INDEL', type=lambda x: x.strip(), default='6,6',
                          help='gap open penalties for deletions and insertions [6,6]',
                          required=False)
align_parser.add_argument('-E', type=lambda x: x.strip(), default='1,1',
                          help='gap extension penalty; a gap of size k cost \'{-O} + {-E}*k\' [1,1]',
                          required=False)
align_parser.add_argument('-L', type=lambda x: x.strip(), default='30,30',
                          help='penalty for 5\'- and 3\'-end clipping [30,30]',
                          required=False)
align_parser.add_argument('-U', type=int, default='17',
                          help='penalty for an unpaired read pair [17]',
                          required=False)
align_parser.add_argument('-p', action='store_true', default=False,
                          help='smart pairing (ignoring in2.fq)',
                          required=False)
align_parser.add_argument('-R', type=str, default=None,
                          help='read group header line such as \'@RG\tID:foo\tSM:bar\' [null]',
                          required=False)
align_parser.add_argument('-H', type=str, default=None,
                          help='insert STR to header if it starts with @; or insert lines in FILE [null]',
                          required=False)
align_parser.add_argument('-j', action='store_true', default=False,
                          help='treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)',
                          required=False)
align_parser.add_argument('-T', type=int, default=80,
                          help='minimum score to output [80], set based on read length',
                          required=False)
align_parser.add_argument('-XA', type=lambda x: x.strip(), default='100,200',
                          help='if there are <INT hits with score >80 percent of the max score, '
                               'output all in XA [5,200]',
                          required=False)
align_parser.add_argument('-M', action='store_true', default=False,
                          help='mark shorter split hits as secondary',
                          required=False)
align_parser.add_argument('-I', type=lambda x: x.strip(), default=None,
                          help='specify the mean, standard deviation (10 percent of the mean if absent), max '
                               '(4 sigma from the mean if absent) and min of the insert size distribution.  '
                               'FR orientation only. [inferred], Float,Float,Int,Int',
                          required=False)


# Add Index Parser Commands

index_parser.add_argument('-G', type=str, required=True,
                          help='Path to reference genome fasta file, fasta file should contain all contigs')
index_parser.add_argument('-DB', type=str, required=True,
                          help='Path to index directory, will create directory if folder does not exist')
index_parser.add_argument('-MR', type=str, default=None, help='Path to bed file of mappable regions. Index will be '
                                                              'built using using masked contig sequence')
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
call_meth_parser.add_argument('-BQ', type=int, default=10, help='Minimum base quality for a base to considered for'
                                                               'methylation calling, default=0')
call_meth_parser.add_argument('-MQ', type=int, default=20, help='Minimum alignment quality for an alignment to be '
                                                                'considered for methylation calling, default=20')
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
sim_parser.add_argument('-CG', default=None, help='Path to CGmap file to generate simulation reference profile')
sim_parser.add_argument('-BR', default=None, help='Path to previously generated BSBolt methylation reference')
sim_parser.add_argument('-MR', type=float, default=0.005, help='Mutation rate')
sim_parser.add_argument('-MI', type=float, default=0.20, help='Mutation indel fraction')
sim_parser.add_argument('-ME', type=float, default=0.20, help='Mutation indel extension probability')
sim_parser.add_argument('-RS', type=int, default=-1, help='Random seed for variant generation')
sim_parser.add_argument('-HA', default=False, action='store_true', help='Haplotype mode, homozygous variants only')
sim_parser.add_argument('-CH', default=True, action='store_false', help='Skip simulation of CH methylation, '
                                                                        'all CH sites unmethylated')
sim_parser.add_argument('-NS', default=True, action='store_false', help='By default observed methylation counts are '
                                                                        'saved, disable this behavior')
sim_parser.add_argument('-SE', type=float, default=0.001, help='Sequencing Error')
sim_parser.add_argument('-NF', type=float, default=0.05, help='Cutoff threshold for amibiguous bases, simulated reads '
                                                              'with a proportion of ambiguous bases above this '
                                                              'threshold will not be output')
sim_parser.add_argument('-FM', type=int, default=400, help='Max fragment size')
sim_parser.add_argument('-IM', type=int, default=50, help='Insert length mean')
sim_parser.add_argument('-SM', type=int, default=50, help='Insert length standard deviation')
sim_parser.add_argument('-verbose', action='store_true', default=False, help='Verbose read simulation')
sim_parser.add_argument('-overwrite', action='store_true', default=False, help='Overwrite previously generated '
                                                                               'simulation database')


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

# Add Sort Parser Args

sort_parser.add_argument('-O', type=str, required=True, help='Sorted bam output')
sort_parser.add_argument('-I', type=str, required=True, help='Input bam')

# add bam indexing args

bam_index.add_argument('-I', type=str, required=True, help='BAM input path')
