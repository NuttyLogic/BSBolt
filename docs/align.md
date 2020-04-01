## Read Alignment


### Adapter Trimming 

It is good practice to trim adapter sequences from sequencing reads before alignment with BSBolt. We recommend using 
[CutAdapt](https://cutadapt.readthedocs.io/en/stable/). 


### BSB Index
Following bisulfite conversion the Watson (sense) and  Crick (anti-sense) DNA strands are no longer complementary. 
Correctly aligning bisulfite converted reads requires the construction of an alignment index containing a bisulfite 
converted sequence for each reference strand; resulting in a larger index. Not all regions of the genome will be 
bisulfite unique, ie reads from these regions will not align uniquely to either the Watson or Crick strand. 
Reads mapping to different reference strands are handled during the alignment step.  

BSB support generation of 3 types of alignment indices. 

1. Whole genome bisulfite sequencing indices 
    - Index is generated for complete Watson and Crick strands
2. Masked alignment indices
    - Sequence outside the reference region is masked in both the Watson and Crick strands before index generation
3. Reduced representation bisulfite sequencing (RRBS) indices
    - The reference sequence is *In silico* restriction enzyme digested to give a list of plausible target regions

**BSB Index Commands**
```shell
  -h, --help            show this help message and exit
  -G                    Path to reference genome fasta file, fasta file should
                        contain all contigs
  -DB                   Path to index directory, will create directory if
                        folder does not exist
  -MR                   Path to bed file of mappable regions. Index will be
                        built using using masked contig sequence
  -rrbs                 Generate a Reduced Representative Bisulfite Sequencing
                        (RRBS) index
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
python3 -m BSBolt Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB 
```

**Masked Alignment Index Generation Example**
```shell
python3 -m BSBolt Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -MR /Tests/TestData/test_wgbs_madking.bed
```

**RRBS Index Generation Example**
```shell
# RRBS Index, MSPI Cut Format, 40bp Lower Fragment Bound, and 400bp Upper Fragment Bound
python3 -m BSBolt Index -G ~/Tests/TestData/BSB_test.fa -DB ~/Tests/TestData/BSB_Test_DB -rrbs -rrbs-cut-format C-CGG -rrbs-lower 40 -rrbs-upper 400
```

### BSB Align
Input reads are modified by first performing *In silico* bisulfite conversion of any unconverted cytosine bases 
(or guanine for the PCR product of bisulfite converted DNA). Converted reads are mapped to a bisulfite alignment 
index and then assessed for bisulfite uniqueness. Alignment mapping quality indicates the uniqueness of the alignment. 
An alignment with several possible alternative mapping locations has a low mapping quality. Additionally, an alignment
 with alternative mappings locations on different strands is assigned a mapping quality of 0, as the alignment is 
 not bisulfite unique. Alignment score thresholds (*-T*) should be set high relative to read length to reduce the frequency of 
observed bisulfite amibiguous reads. Valid alignments are further modified so all Watson reads are reported as sense (+) 
reads and all Crick reads are reported as anti-sense (-) reads.  

**BSB Align Commands**
```shell
  -F1           Path to fastq 1
  -F2           Path to fastq 2
  -UN           Library Undirectional, Consider PCR products of bisulfite converted DNA
  -O            Path to Output Prefix
  -G            Path to BSBolt Database
  -t            Number of bwa threads [1]
  -k            Minimum seed length [19]
  -w            Band width for banded alignment [100]
  -d            off-diagonal X-dropoff [100]
  -r            look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
  -y            seed occurrence for the 3rd round seeding [20]
  -c            skip seeds with more than INT occurrences [500]
  -D            drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
  -W            discard a chain if seeded bases shorter than INT [0]
  -m            perform at most INT rounds of mate rescues for each read [50]
  -S            skip mate rescue
  -P            skip pairing; mate rescue performed unless -S also in use
  -A            score for a sequence match, which scales options -TdBOELU unless overridden [1]
  -B            penalty for a mismatch [4]
  -INDEL        gap open penalties for deletions and insertions [6,6]
  -E            gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
  -L            penalty for 5'- and 3'-end clipping [30,30]
  -U            penalty for an unpaired read pair [17]
  -p            smart pairing (ignoring in2.fq)
  -R            read group header line such as '@RG ID:foo SM:bar' [null]
  -H            insert STR to header if it starts with @; or insert lines in FILE [null]
  -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)
  -T            minimum score to output [80], set based on read length
  -XA           if there are <INT hits with score >80 percent of the max score, output all in XA
                [5,200]
  -M            mark shorter split hits as secondary
  -I            specify the mean, standard deviation (10 percent of the mean if absent), max (4
                sigma from the mean if absent) and min of the insert size distribution. FR
                orientation only. [inferred], Float,Float,Int,Int
```
**Paired End Alignment**
```shell
# Paired End Alignment Using Default Commands
python3 -m BSBolt Align -DB ~/Tests/TestData/BSB_Test_DB -F1 ~/Tests/TestSimulations/BSB_pe_meth_1.fastq -F2 ~/Tests/TestSimulations/BSB_pe_meth_2.fastq -O ~/Tests/BSB_pe_test 
```

**Single End Alignment**
```shell
# Single End Alignment Using Default Commands
python3 -m BSBolt Align -DB ~/Tests/TestData/BSB_Test_DB -F1 ~/Tests/TestSimulations/BSB_pe_meth_1.fastq -O ~/Tests/BSB_pe_test 
```

**Alignment Output**
Aligned reads are output as a BAM file. BSBolt formats the aligned reads to convey bisulfite specific alignment information. 
Watson reads alignments are reported on the positive strand and Crick alignments on the negative strand. This information is conveyed 
in the sam flag, and in a BSBolt sam tag. All tags are listed in the table below.

| Tag | Type    | Description | Example |
| :---: | :---:    | :---: | :---: |
|YS    |Z (string)| Mapping strand (C=Crick, W=Watson) and alignment conversion pattern (C2T or G2A) | YS:Z:W_C2T (ie. Watson_Cytosine.to.Thymine)|  
|YC    |i (integer)| Bisulfite ambiguous alignments, 1 if ambiguous, not reported if non-ambiguous| YC:i:1 