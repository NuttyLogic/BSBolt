# BSBolt (BiSulfite Bolt)
## A fast and safe read alignment platform for bisulfite sequencing data
*BSBolt is in pre-release, all documentation and feature usage is not finalized*

BSBolt is fast and safe bisulfite sequencing processing platform. BSBolt offers support for bisulfite sequencing 
read simulation, read alignment, methylation calling, and methylation matrix assembly. BSBolt is 
an evolution of [BSSeeker2](https://github.com/BSSeeker/BSseeker2); re-imagined for effective processing of large 
bisulfite sequencing experiments.   

## BSBolt Content
- [BSBolt (BiSulfite Bolt)](#bsbolt-bisulfite-bolt)
  - [A fast and safe read alignment platform for bisulfite sequencing data](#a-fast-and-safe-read-alignment-platform-for-bisulfite-sequencing-data)
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
Representative Bisulfite Sequencing* (RRBS), and targeted bisulfite sequencing data. BSBolt utililzies 
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) as its alignment engine and 
[ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) for simulation of 
Illumina reads. 

### Installation
BSBolt requries python3.6 or greater to run, with the external *pysam*, *numpy*, and *tqdm* packages. BSBolt is 
distributed with [Bowtie2]( http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
and [ART Read Simulation Tools](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) 
under the terms of the [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.html). BSBolt is implemented for Unix like
operating systems (MacOS, Linux); Windows support is untested. 

## Read Alignment

### Data Preprocessing
Adapter sequences should be trimmed before alignment with BSBolt. We recommend using 
[CutAdapt](https://cutadapt.readthedocs.io/en/stable/) to trim sequencing adapters. If aligning un-trimmed FASTQ files 
align reads using Bowtie2 local mode for best performance.
### BSB Index 
Read alignment requires a BSBolt index to proceed. Alignment indexes for *RRBS* and 
*WGBS* must be generated separately.
 
**BSB Index Commands**
```bash
-h, --help          show this help message and exit
-G                  Path for reference genome fasta file, fasta file
                    should contain all contigs
-DB                 Path to index directory, will create directory if
                    folder does not exist
-BT2                Path to bowtie2 executable
-BT2-p              Number of threads for Bowtie2 to use
-rrbs               Generate Reduced Representative Bisulfite Sequencing Index
-rrbs-cut-format    RRBS_CUT_FORMAT
                    Cut format to use for generation of RRBS database,
                    default= C-CGG (MSPI), input multiple enzymes as a
                    comma seperate string, C-CGG,C-CGG,...
-rrbs-lower         Lower bound fragment size to consider RRBS
                    indexgeneration, default = 40
-rrbs-u             Upper bound fragment size to consider RRBS
                    indexgeneration, default = 500
```
**WGBS Index Generation Example**
```bash
# WGBS Index with 4 BT2 Threads
python3 BSBolt.py Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -BT2 bowtie2 -BT2-p 4
```

**RRBS Index Generation Example**
```bash
# RRBS Index Using 4 BT2 Threads, MSPI Cut Format, 40bp Lower Fragment Bound, and 400bp Upper Fragment Bound
python3 BSBolt.py Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -BT2 bowtie2 -BT2-p 4 -rrbs -rrbs-cut-format C-CGG -rrbs-lower 40 -rrbs-upper 400
```


### BSB Align
Following index generation alignment of *RRBS* and *WGBS* follows the same pipeline.

**BSB Align Commands**
```bash
-h, --help            show this help message and exit
-F1                   Path to fastq 1
-F2                   Path to fastq 2
-U                    Library undirectioinal, default=True
-BT2                  Path to bowtie2 aligner
-O                    Path to Output Prefix
-DB                   Path to BSBolt Database
-CP                   Proportion threshold to label read not fully converted
-CT                   Number of mCH that must be observed to label a read
                      unconverted
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
```bash
# Paired End Alignment Using Default Commands
python3 BSBolt.py Align -DB ~/Tests/TestData/BSB_Test_DB -BT2 bowtie2 -F1 ~/Tests/TestSimulations/BSB_pe_meth_1.fastq -F2 ~/Tests/TestSimulations/BSB_pe_meth_2.fastq -O ~/Tests/BSB_pe_test -S
```

**Single End Alignment**
```bash
# Single End Alignment Using Default Commands
python3 BSBolt.py Align -DB ~/Tests/TestData/BSB_Test_DB -BT2 bowtie2 -F1 ~/Tests/TestSimulations/BSB_pe_meth_1.fastq -O ~/Tests/BSB_pe_test -S
```

## BSB CallMethylation
To ensure efficient methylation calling BSB-Call-Methylation only supports sorted BAM files as the alignment input.
**BSB CallMethylation Commands**
```bash
-h, --help       show this help message and exit
-I               Input BAM, input file must be in BAM format
-DB              Path to index directory
-O               Output prefix
-remove-ccgg     Remove methylation calls in ccgg sites, default=False
-verbose         Verbose Output, default=False
-text            Output plain text files, default=False
-remove-sx       Remove methylation calls from reads marked as incompletely
                 by BSBolt Align, default=True
-ignore-overlap  Only consider higher quality base when paired end reads
                 overlap, default=True
-max             Max read depth to call methylation
-min             Minimum read depth required to report methylation site
-t               Number of threads to use when calling methylation values
-min-qual        Minimum base quality for a base to considered for methylation calling
```

**Methylation Calling**
```bash
# Methylation Calling with 2 threads, 
python3 BSBolt.py CallMethylation -I ~/Tests/BSB_pe_test.sorted.bam -O ~/Tests/BSB_pe_test -DB ~/Tests/TestData/BSB_Test_DB -t 2 -verbose > methylation_stats.txt
```
**Output Files**
BSB CallMethylation outputs ATCGmap, CGmap, and wig files by default.

## BSB Simulate
**BSB Simulate Commands**
```bash
-h, --help  show this help message and exit
-G          Path for reference genome fasta file, fasta file should contain
            all contigs
-A          Path to ART executable
-O          Output prefix
-PE         Simulate Paired End Reads, default Single End
-RL         Simulated Read Lenghth
-RD         Simulated Read Depth
-U          Simulate Undirectional Reads, default=Directional
```
**Simulate Paired End, Undirectional Methylation Reads**
```bash
python3 BSBolt.py Simulate -G ~/Tests/TestData/BSB_test.fa -A ~/art_bin_MountRainier/art_illumina -O ~/Tests/TestSimulations/BSB_pe -U -PE
```

## Methylation Matrix Assembly
**BSBolt AggregateMatrix Commands**
```bash
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
