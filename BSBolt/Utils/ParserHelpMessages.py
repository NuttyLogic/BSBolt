alignment_help = '''
BSBolt Align -F1 {fastq1} -DB {BSBolt db} -O {output}

-h, --help    show this help message and exit

Input / Output Options:
  -F1 File     path to fastq 1
  -F2 File     path to fastq 2 [null]
  -O File      output Prefix
  -DB File     path to BSBolt database
  -R Str       read group header line such as '@RG ID:foo SM:bar' [null]
  -H Str       insert STR to header if it starts with @; or insert lines in FILE [null]
  -XA Int,Int  if there are <INT hits with score >80 percent of the max score, output all in XA [100,200]
  -DR Float    drop ratio for alternative hits reported in XA tag, 
               for best bisulfite alignment performance set at or above default [0.95]
  -p           smart pairing (ignoring in2.fq)
Scoring Options
  -A Int       score for a sequence match, which scales options -TdBOELU unless overridden [1]
  -B Int       penalty for a mismatch [4]
  -INDEL Int   gap open penalties for deletions and insertions [6,6]
  -E Int       gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
  -L Int,Int   penalty for 5'- and 3'-end clipping [30,30]
  -U Int       penalty for an unpaired read pair [17]
Bisulfite Options
  -UN          library undirectional, ie. consider PCR products of bisulfite converted DNA
  -CP Float    CH conversion proportion threshold [0.5]
  -CT Int      number of CH sites needed to assess read conversion
  -SP Float    substitution threshold for read bisulfite conversion patterns (ie C2T, G2A) [0.1]
               for undirectional libraries the substitution pattern with the fewer number of 
               substitutions relative to the total read length (if < threshold) is aligned preferentially
Algorithm Options
  -t Int       number of bwa threads [1]
  -k Int       minimum seed length [19]
  -w Int       band width for banded alignment [100]
  -d Int       off-diagonal X drop off [100]
  -r Float     look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
  -y Int       seed occurrence for the 3rd round seeding [20]
  -c Int       skip seeds with more than INT occurrences [500]
  -D Float     drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
  -W Int       discard a chain if seeded bases shorter than INT [0]
  -m Int       perform at most INT rounds of mate rescues for each read [50]
  -S           skip mate rescue
  -P           skip pairing; mate rescue performed unless -S also in use
  -j           ignore ALT contigs
  -T Int       minimum score to output [80], set based on read length
  -M           mark shorter split hits as secondary
  -I Fl,Fl,Int,Int    specify the mean, standard deviation (10 percent of the mean if absent), 
                      max (4 sigma from the mean if absent) and min of the insert size distribution. 
                      FR orientation only. [inferred]

'''

index_help = '''
BSBolt Index -G {fasta reference} -DB {database output}

-h, --help     show this help message and exit

Index Options:
  -G File      path to reference genome fasta file, fasta file should contain all contigs
  -DB File     path to index directory, alignment index generated inside existing directory 
               or new directory made if directory doesn't exist
  -B Int       block size for bwtsw algorithm, 
               increasing block will speed up indexing and increase memory consumption [10000000]
  -MR File     path to bed file of mappable regions. 
               Index will be built using using masked contig sequence [null]
  -IA          ignore alt contigs when constructing alignment index
  -rrbs        generate a Reduced Representative Bisulfite Sequencing (RRBS) index
  -rrbs-cut-format Str    Cut format to use for generation of RRBS database, [C-CGG] MSPI, 
                          input multiple enzymes as a comma separate string, C-CGG,C-CGG,...
  -rrbs-lower Int      lower bound fragment size to consider RRBS index generation [40]
  -rrbs-upper Int      upper bound fragment size to consider RRBS index generation [500]
'''

meth_help = '''
BSBolt Module CallMethylation -I {input.bam} -DB {BSBolt DB} -O {output prefix}

-h, --help     show this help message and exit

Input / Output Options:
  -I File         input BAM, input file must be in BAM format with index file
  -DB File        path to index directory
  -O File         output prefix
  -text           output plain text files [False]
  -CG             only output CpG sites in CGmap file [False]
Algorithm Options:
  -remove-ccgg    remove methylation calls in ccgg sites [False]
  -verbose        verbose Output [False]
  -ignore-ov      only consider higher quality base when paired end reads overlap [True]
  -max Int        max read depth to call methylation [8000]
  -min Int        minimum read depth required to report methylation site [10]
  -t Int          methylation calling threads [1]
  -BQ Int         minimum base quality [10]
  -MQ Int         minimum alignment quality [20]
  -IO             ignore orphans reads, (not proper read pair)
'''

aggregate_help = '''
BSBolt AggregateMatrix -F {file1.CGmap,file2.CGmap,...} -O {output_matrix.txt}

-h, --help              show this help message and exit

Options:
  -F File,File,.        comma separated list of CGmap files, 
                        or path to text file with list of line separated CGmap files
  -S Str,Str,.          comma separated list of samples labels. If sample labels are not provided sample labels 
                        are extracted from CGmap files. Can also pass path to txt for line separated sample labels.
  -min-coverage Int     minimum site read depth coverage 
  -min-sample Float     proportion of samples that must have a valid site (above minimum coverage threshold)
  -O File               Aggregate matrix output path
  -CG                   Only output CG sites
  -verbose              Verbose aggregation
  -t Int                Number of threads to use when assembling matrix
  -count                Output a count matrix with count of methylated cytosines and total observed cytosines
'''

sim_help = '''
BSBolt Simulate -G {genome.fa} -O {output_directory}

-h, --help  show this help message and exit

Input / Output Options:
  -G File     path for reference genome fasta file
  -O File     output prefix
  -CG File    path to CGmap file reference profile [Null]
  -overwrite  overwrite previously generated simulation database
  -BR File    Path to previously generated BSBolt methylation reference (directory)
  -NS         don't output simulated methylation counts
  -verbose    verbose read simulation
Algorithm Options:
  -PE         simulate Paired End Reads, default Single End
  -RL Int     simulated Read Length [125]
  -RD Int     simulated Read Depth [20]
  -U          simulate undirectional reads, (bisulfite converted reference strands and PCR products)
  -MR Float   mutation rate [0.005]
  -MI Float   mutation indel fraction [0.20]
  -ME Float   mutation indel extension probability [0.20]
  -RS Int     random seed for variant generation [-1]
  -HA         haplotype mode, homozygous variants only
  -CH         skip simulation of CH methylation, all CH sites unmethylated
  -SE Float   sequencing Error [0.001]
  -NF Float   cutoff threshold for ambiguous bases, simulated reads with a proportion of ambiguous 
              bases above this threshold will not be output [0.05]
  -FM Int     max fragment size [400]
  -IM Int     insert length mean [50]
  -SM Int     insert length standard deviation [50]
'''

impute_help = '''
BSBolt Impute -M {BSBolt_matrix.txt} -O {imputed_matrix.txt}

-h, --help  show this help message and exit

Options:
  -M File     path to BSBolt matrix file
  -B Int      imputation sample batch size kNN imputation, by default the all of the samples 
              will be processed as a single batch
  -W Int      sliding window size for imputation [3000000]
  -k Int      number of neighbors to use for imputation [5]
  -t Int      number of threads available for imputation [1]
  -verbose    verbose imputation
  -O File     output path for imputed matrix
  -R          randomize batches
'''