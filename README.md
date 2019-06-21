# BSBolt (BiSulfite Bolt)
## A fast and safe read alignment platform for bisulfite sequencing data

BiSuflite Bolt (BSBolt); a fast and scalable bisulfite sequencing analysis platform. BSBolt is an integrated analysis 
platform that offers support for bisulfite sequencing read simulation, alignment, methylation calling, data aggregation, 
and data imputation. BSBolt has been validated to work with a wide array of bisulfite sequencing data, 
including whole genome bisulfite sequencing (WGBS), reduced representative bisulfite sequencing data (RRBS), 
and targeted methylation sequencing data. 
  
## Documentation

The complete BSBolt documentation can be found at, [bsbolt.readthedocs.io](bsbolt.readthedocs.io).

## BSBolt Content
- [BSBolt (BiSulfite Bolt)](#bsbolt-bisulfite-bolt)
  - [BSBolt Content](#bsbolt-content)
  - [Overview](#overview)
    - [Installation](#installation)
  - [Read Alignment](#read-alignment)
    - [Data Preprocessing](#data-preprocessing)
    - [BSB-Index](#bsb-index)
    - [BSB-Align](#bsb-align)
  - [BSB-Call-Methylation](#bsb-call-methylation)
  - [BSB-Simulate](#bsb-simulate)
  - [Methylation Matrix Assembly](#methylation-matrix-assembly)

## Overview
BSBolt provides support for processing and simulation of *Whole Genome Bisulfite Sequencing* (WGBS), *Reduced 
Representative Bisulfite Sequencing* (RRBS), and targeted bisulfite sequencing data. BSBolt utilizes
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) as its alignment engine and 
[ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) for simulation of 
Illumina reads. 

### Installation
BSBolt requires python3.6 or greater to run, with the external *pysam*, *numpy*, and *tqdm* packages. BSBolt is 
distributed with [Bowtie2]( http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
and [ART Read Simulation Tools](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) 
under the terms of the [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.html). BSBolt is implemented for Unix like
operating systems (MacOS, Linux); Windows support is untested. 

**Calling BSBolt Modules**
Invoke BSBolt modules by calling *BSBolt.py* and selecting the appropriate module. Module help command can be invoked 
by calling *-h* with each module. 
```
python3 BSBolt.py Module

Align               Alignment Module
Index               Index Generation Module
CallMethylation     Methylation Calling Module
AggregateMatrix     CGmap Matrix Aggregation Module
Simulate            BSBolt Illumina Read Simulation Module
Impute              kNN Imputation Module
```

## Read Alignment

### Data Pre-Processing
Adapter sequences should be trimmed from sequencing reads before alignment with BSBolt. We recommend using 
[CutAdapt](https://cutadapt.readthedocs.io/en/stable/) to trim sequencing adapters. If aligning un-trimmed FASTQ files,
alignment can be performed using Bowtie2 local mode to increase mappability; however, due the low complexity of bisulfite
sequencing libraries this is not recommended. 
### BSB Index 
Alignment with BSBolt requires generation of an alignment index. BSBbolt Index generates multiple Bowtie2 alignment indices in 
addition to processed reference files for downstream analysis. 

BSB support generation of 3 types of alignment indices. 
1. Whole genome bisulfite sequencing indices 
2. Masked alignment indices for targeted bisulfite sequencing 
3. Reduced representation bisulfite sequencing indices (*In silico* restriction enzyme digested)
 
**BSB Index Commands**
```
  -h, --help            show this help message and exit
  -G                    Path for reference genome fasta file, fasta file
                        should contain all contigs
  -DB                   Path to index directory, will create directory if
                        folder does not exist
  -MR                   Path to bed file of mappable regions, build Bowtie2
                        reference using masked contig sequence
  -BT2                  Path to bowtie2 executable, default = bundled bowtie2
  -BT2-p                Number of threads for Bowtie2 to use
  -rrbs                 Generate Reduced Representative Bisulfite Sequencing
                        Index
  -rrbs-cut-format      Cut format to use for generation of RRBS database,
                        default= C-CGG (MSPI), input multiple enzymes as a
                        comma seperate string, C-CGG,C-CGG,...
  -rrbs-lower           Lower bound fragment size to consider RRBS
                        indexgeneration, default = 40
  -rrbs-upper           Upper bound fragment size to consider RRBS
                        indexgeneration, default = 500
```
**WGBS Index Generation Example**
```
# WGBS Index with 4 BT2 Threads
python3 BSBolt.py Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -BT2-p 4
```

**Masked Alignment Index Generation Example**
```
# WGBS Index with 4 BT2 Threads
python3 BSBolt.py Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -BT2-p 4 -MR /Tests/TestData/test_wgbs_madking.bed
```

**RRBS Index Generation Example**
```
# RRBS Index Using 4 BT2 Threads, MSPI Cut Format, 40bp Lower Fragment Bound, and 400bp Upper Fragment Bound
python3 BSBolt.py Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -BT2-p 4 -rrbs -rrbs-cut-format C-CGG -rrbs-lower 40 -rrbs-upper 400
```


### BSB Align
Read alignment is performed using the BSBolt Align module. During read alignment, sequenced bases that were not converted 
during bisulfite treatment are converted *in silico* and aligned to a bisulfite converted Watson (sense) and a Crick (anti-sense) 
reference genome. Aligned reads are selected based mapping criteria provided by the user and incorporated into a final alignment file.

**BSB Align Commands**
```
-h, --help            show this help message and exit
-F1                   Path to fastq 1
-F2                   Path to fastq 2
-NC                   Aligned unconverted bisulfite reads
-OU                   Output unmapped reads
-U                    Library undirectioinal, default=True
-BT2                  Path to bowtie2 aligner
-O                    Path to Output Prefix
-DB                   Path to BSBolt Database
-M                    Read mismatch threshold, reads with mismatches greater
                      than threshold will be discarded
-S                    Position Sort Output Bam, default=False
-BT2-local            Bowtie2 local alignment or end-to-end, default end-to-
                      end
-BT2-D                Bowtie2 number of consecutive seed extension attempts
                      that can fail before Bowtie2 move on
-BT2-k                Bowtie2 alignment search limit
-BT2-p                Number of threads for Bowtie2 to use
-BT2-L                Length of subseeds during alignment
-BT2-score-min        Bowtie2 scoring function
-BT2-I                Bowtie2, minimum fragment length for a valid paired-
                      end alignment
-BT2-X                Bowtie2, maximum fragment length for a valid paired-
                      end alignment
```
**Paired End Alignment**
```
# Paired End Alignment Using Default Commands
python3 BSBolt.py Align -DB ~/Tests/TestData/BSB_Test_DB -F1 ~/Tests/TestSimulations/BSB_pe_meth_1.fastq -F2 ~/Tests/TestSimulations/BSB_pe_meth_2.fastq -O ~/Tests/BSB_pe_test -S
```

**Single End Alignment**
```
# Single End Alignment Using Default Commands
python3 BSBolt.py Align -DB ~/Tests/TestData/BSB_Test_DB -F1 ~/Tests/TestSimulations/BSB_pe_meth_1.fastq -O ~/Tests/BSB_pe_test -S
```

## Methylation Calling 

Methylation values, the proportion of methylated bases to the total number of observed bases, are called using the BSBolt CallMethylation module. 
Methylation values for guanine nucleotides are made for reads aligned to the Crick (anti-sense) strands and calls for cytosine nucleotides are made for 
reads aligned to the Watson (sense) strand. Methylation calls for WGBS and targeted bisuflite sequencing can be improved by removing PCR duplicate reads 
before calling methylation. It is not recommend to remove duplicated for RRBS as the sequencing reads will often share the same mapping coordinates 
due to enzymatic digestion. 

### Methylation Calling Pre-Processing

We recommend using [samtools](http://www.htslib.org/) to remove duplicates. Due to the structure of the integrated alignment file the *samtools fixmate* option must be disabled. 
Removing PCR duplicates using paired sequencing reads will give better results; removal using single end reads can be overly aggressive and should be performed on a case by case basis. 

```
# fixmates to prepare for duplicate removal, use -p to disable proper pair check
samtools samtools fixmate -p BSB_pe_test.bam BSB_pe_test.fixmates.bam 
# remove duplicate reads
samtools markdup -r BSB_pe_test.fixmates.bam BSB_pe_test.dup.bam
# sort bam by coordinates for methylation calling
samtools sort BSB_pe_test.dup.bam BSB_pe_test.sorted.bam
```

### BSBolt CallMethylation

The running time of the methylation calling module can be greatly reduced by only calling methylation for CG sites. By default compressed CGmap files are output from methylation calling. 
ATCGmap files can also be output if downstream analysis requires, but this is disabled by default

**BSB CallMethylation Commands**
```
  -h, --help          show this help message and exit
  -I                  Input BAM, input file must be in BAM format
  -DB                 Path to index directory
  -O                  Output prefix
  -remove-ccgg        Remove methylation calls in ccgg sites,default=False
  -verbose            Verbose Output, default=False
  -text               Output plain text files, default=False
  -ignore-overlap     Only consider higher quality base when paired end reads
                      overlap, default=False
  -max                Max read depth to call methylation
  -min                Minimum read depth required to report methylation site
  -t                  Number of threads to use when calling methylation values
  -min-qual           Minimum base quality for a base to considered
                      for methylation calling, default=0
  -CG                 Only output CpG sites in CGmap file
  -ATCG               Output ATCGmap file
```

**Methylation Calling**
Methylation calling is performed by counting the number of bisulfite converted bases relative to the number of reads 
observed at each cytonsine. Relative to the reference genome methylation status at a cytosine and guanine 
can only be called using reads mapped to Watson and Crick strands respectively. 
```
# Methylation Calling with 2 threads, 
python3 BSBolt.py CallMethylation -I ~/Tests/BSB_pe_test.sorted.bam -O ~/Tests/BSB_pe_test -DB ~/Tests/TestData/BSB_Test_DB -t 2 -verbose > methylation_stats.txt
```
**Output Files**
CGmap is a tab separated txt format desribing the methylation status of a cytosine. 
1. Chromosome
2. Nucleotide, C for reads mapped to the Watson (sense) strand and G for reads mapped to the Crick (anti-sense) strand
3. Position, base-pairs from start
4. Context, three base pair methylation context
5. Sub-Context, two base pair methylation context
5. Methylation Value, proportion of methylation reads to total reads
6. Methylation Bases, methylated nucleotides observed
7. All Bases, total number of nucleotides observed at the mapping position

```text
chrom   nucleotide  position   context sub-context  methylation_value methylated_bases  all_bases 
chr11	G	422436	CHH	CC	0.091	1	10
chr12	G	389290	CHH	CT	0.0	0	10
chr13	G	200552	CHH	CT	0.0	0	10
chr11	C	142826	CG	CG	0.0	0	10
```

## Methylation Matrix Assembly

BSBolt AggregrateMatrix takes a list of CGmap, compressed or uncompressed, and assembles a consensus methylation matrix. Methylated sites that 
pass a read depth threshold and are present in a set proportion of samples are included in the matrix. 

**BSBolt AggregateMatrix Commands**
```
optional arguments:
  -h, --help            show this help message and exit
  -F                    Comma separated list of CGmap file paths, or path to
                        text file with list of line separated CGmap file paths
  -S                    Comma separated list of samples labels. If sample
                        labels are not provided sample labels are extracted
                        from CGmap file paths. Can also pass path to txt for
                        line separated sample labels.
  -min-coverage         Minimum site read depth coverage for a site to be
                        included in the aggregate matrix
  -min-sample           Proportion of samples that must have a valid site
                        (above minimum coverage threshold), for a site to
                        beincluded in the aggregate matrix.
  -O                    Aggregate matrix output path
  -CG                   Only output CG sites
  -verbose              Verbose aggregation
```
**Aggregate Matrix Default Settings**

```bash
python3 BSBolt.py AggregateMatrix -F cgmap_1,cgmap_2,cgmap_3 -O ~/test_matrix.txt
```
**Aggregate Matrix Default Settings - File List**

```bash
python3 BSBolt.py AggregateMatrix -F cgmap_file_list.txt -O ~/test_matrix.txt
```

**Aggregate Matrix Default Settings - File List, Sample Labels, Verbose**

```bash
python3 BSBolt.py AggregateMatrix -F cgmap_file_list.txt -S sample1,sample2,sample3 -O ~/test_matrix.txt -verbose
```

## Methylation Value Imputation

BSBolt Impute leverages the correlation structure between neighboring CpG sites to impute missing values through the use of a kNN sliding window.  Within each window the nearest neighbors are calculated using Euclidean distance, 
and the average value of k nearest neighbors is used to impute the null methylation value. To efficiently scale the algorithm, imputation can be performed in batches. 

```
  -h, --help  show this help message and exit
  -M          Path to BSB matrix file
  -B          Imputation sample batch size kNN imputation, by default the all
              of the samples will be processed as a single batch
  -W          Sliding window size for imputation
  -k          Number of neighbors to use for imputation, default = 5
  -t          Number of threads available for imputation
  -verbose    Verbose output
  -O          Output path for imputed matrix
  -R          Randomize batches
```  

**Impute No Batches**
```bash
python3 BSBolt.py Impute -M ~/test_matrix.txt -W 100000 -k 3 -t 4 -O ~/test_matrix.impute.txt
```

**Batch Imputation**
```bash
python3 BSBolt.py Impute -M ~/test_matrix.txt -W 100000 -k 3 -t 4 -O ~/test_matrix.impute.txt -B 10 -R
```

## Bisulfite Sequencing Simulation

Simulating bisuflite reads is performed by assiging each cytosine and guanine in the user provided reference file a 
methylation value. This done based on an input methylation reference file (.CGmap or BSBolt Reference Directory) or 
sampling from a binomial distribution for CpG and CH sites. After a methylation value is assigned for cytosine and 
guanine individual reads simulated using ART are modified according to the methylation value of the nucleotides in the 
read. The methylation value represent the probability a base with be methylated, and as a result there random noise is 
introduced to more accurately represent bisulfite sequencing reads. 

**BSB Simulate Commands**
```
  -h, --help  show this help message and exit
  -G          Path for reference genome fasta file, fasta file should contain
              all contigs
  -A          Path to ART executable, default = bundled ART
  -O          Output prefix
  -PE         Simulate Paired End Reads, default Single End
  -RL         Simulated Read Lenghth
  -RD         Simulated Read Depth
  -U          Simulate Undirectional Reads, default=Directional
  -RC         Path to CGmap file to generate simulation reference profile
  -RO         Methylation reference output directory, default = output path
  -BR         Path to previously generate BSB simulation reference
  -IR1        Read 1 insertion rate
  -IR2        Read 2 insertion rate
  -DR1        Read 1 deletion rate
  -DR2        Read 2 deletion rate
  -NF         Cutoff theshold for a read with gaps, -, or ambiguous bases, N.
              Reads below the threshold will not be output
  -M          Mean paired end fragment size
  -SM         Paired end fragment length distribution standard deviation
  -SS         Sequencing system, default HS25
              ART Sequencing System Codes
              * HS10 - HiSeq 1000 (100bp) 
              * HS20 - HiSeq 2000 (100bp) 
              * HS25 - HiSeq 2500 (125bp, 150bp) 
              * HSXn - HiSeqX PCR free (150bp) 
              * HSXt - HiSeqX TruSeq (150bp) 
              * MinS - MiniSeq TruSeq (50bp) 
              * MSv1 - MiSeq v1 (250bp)
              * MSv3 - MiSeq v3 (250bp) 
              * NS50 - NextSeq500 v2 (75bp)
  -Q1         Optional read 1 quality profile, generated using ART Illumina
              read profiler
  -Q2         Path to read 2 quality profile

```
**Simulate Paired End, Undirectional Methylation Reads**
```bash
python3 BSBolt.py Simulate -G ~/Tests/TestData/BSB_test.fa -O ~/Tests/TestSimulations/BSB_pe -U -PE
```