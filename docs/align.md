## Read Alignment
BSBolt support single and paired end read alignment. Reads can be aligned in end-to-end, where the entire reads is 
aligned without read trimming, or local mode where a read can be trimmed during alignment. Due to the low complexity of 
bisulfite sequencing data end-to-end alignment is generally preferred. If aligning un-trimmed FASTQ files,
alignment can be performed using Bowtie2 local mode; however, due the low complexity of bisulfite
sequencing libraries this is not recommended.  


### Adapter Trimming 

Adapter sequences should be trimmed from sequencing reads before alignment with BSBolt. We recommend using 
[CutAdapt](https://cutadapt.readthedocs.io/en/stable/) to trim sequencing adapters. Illumina adapter sequences commonly 
used for bisulfite sequencing are listed below. 

**Standard Illumina TruSeq Adapter Sequences**

```test
Read 1
    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

Read 2
    AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```

### BSB Index
Following bisulfite conversion the Watson (sense) and  Crick (anti-sense) DNA strands are no longer complementary. 
Correctly aligning bisulfite converted reads requires the construction of an alignment index containing bisulfite 
converted version of each reference strand. The combined reference sequence is twice as large and results in a larger 
alignment index. Reads mapping to different reference strands are handled during the alignment step.  

BSB support generation of 3 types of alignment indices. Local alignment is discourage for when aligning reads against 
as masked alignment index or RRBS index. 

1. Whole genome bisulfite sequencing indices 
    - Index is generated for complete Watson and Crick strands
2. Masked alignment indices for targeted bisulfite sequencing 
    - Sequence outside the reference region is masked in both the Watson and Crick strands before index generation
3. Reduced representation bisulfite sequencing (RRBS) indices
    - The reference sequence is *In silico* restriction enzyme digested to give a list of plausible target region

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
Read alignment is performed by modifying the input sequence by converting any bases not converted during bisulfite treatment. Depending on the 
read bases are converted from cytosine to thymine, or guanine to adenine for the complement of bisulfite converted DNA. Reads from DNA that 
has not been fully bisulfite converted will generally map to both reference strands. In this case the reads are not considered bisulfite unique 
and reported as unmapped reads. Reads produced from the complement of bisulfite converted DNA can also be mapped. In this case the 
pattern of *in silico* base substitution is reversed. Read mapping is handled by Bowtie2. Bowtie2 setting may be changed by passing the appropriate 
argument through the BSBolt Align module. 

A mapped read is considered a valid read alignment if the number of mismatches between the reference and read is below the user set threshold, *-M*, and the 
read is bisulfite unique. Bisulfite unique reads only map to one reference strand implying the read is fully bisulfite converted. Valid reads are further modified 
so all Watson reads are reported as sense (+) reads and all Crick reads are reported as anti-sense (-) reads.  


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
-BT2-local            Bowtie2 local alignment, default end-to-end
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