## Read Alignment

### Data Pre-Processing
Adapter sequences should be trimmed from sequencing reads before alignment with BSBolt. We recommend using 
[CutAdapt](https://cutadapt.readthedocs.io/en/stable/) to trim sequencing adapters. If aligning un-trimmed FASTQ files,
alignment can be performed using Bowtie2 local mode to increase mappability; however, due the low complexity of bisulfite
sequencing libraries this is not recommended. 
### BSB Index 
Alignment with BSBolt requires generation of an alignment index. BSBbolt Index generates a bisulfite converted reference sequence. 

BSB support generation of 3 types of alignment indices. 

1. Whole genome bisulfite sequencing indices 
2. Masked alignment indices for targeted bisulfite sequencing 
3. Reduced representation bisulfite sequencing indices (*In silico* restriction enzyme digested)
 
**BSB Index Commands**
```shell
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
```shell
# WGBS Index with 4 BT2 Threads
python3 BSBolt.py Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -BT2-p 4
```

**Masked Alignment Index Generation Example**
```shell
# WGBS Index with 4 BT2 Threads
python3 BSBolt.py Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -BT2-p 4 -MR /Tests/TestData/test_wgbs_madking.bed
```

**RRBS Index Generation Example**
```shell
# RRBS Index Using 4 BT2 Threads, MSPI Cut Format, 40bp Lower Fragment Bound, and 400bp Upper Fragment Bound
python3 BSBolt.py Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -BT2-p 4 -rrbs -rrbs-cut-format C-CGG -rrbs-lower 40 -rrbs-upper 400
```

### BSB Align
Read alignment is performed using the BSBolt Align module. During read alignment, sequenced bases that were not converted 
during bisulfite treatment are converted *in silico* and aligned to a bisulfite converted Watson (sense) and a Crick (anti-sense) 
reference genome. Aligned reads are selected based mapping criteria provided by the user and incorporated into a final alignment file.

**BSB Align Commands**
```shell
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
```shell
# Paired End Alignment Using Default Commands
python3 BSBolt.py Align -DB ~/Tests/TestData/BSB_Test_DB -F1 ~/Tests/TestSimulations/BSB_pe_meth_1.fastq -F2 ~/Tests/TestSimulations/BSB_pe_meth_2.fastq -O ~/Tests/BSB_pe_test -S
```

**Single End Alignment**
```shell
# Single End Alignment Using Default Commands
python3 BSBolt.py Align -DB ~/Tests/TestData/BSB_Test_DB -F1 ~/Tests/TestSimulations/BSB_pe_meth_1.fastq -O ~/Tests/BSB_pe_test -S
```