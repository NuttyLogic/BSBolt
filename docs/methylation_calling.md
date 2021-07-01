Methylation values, the proportion of methylated bases to the total number of observed bases, are called using the
BSBolt CallMethylation module. Methylation values for guanine nucleotides are made for reads aligned to the Crick
(anti-sense) strands and calls for cytosine nucleotides are made for reads aligned to the Watson (sense) strand.
Methylation calls for WGBS and targeted bisulfite sequencing can be improved by removing PCR duplicate reads before
calling methylation. Marking duplicate alignments with RRBS data is not recommended as the sequencing reads will often share
the same mapping coordinates due to enzymatic digestion.

### Methylation Calling Pre-Processing

We recommend using [samtools](http://www.htslib.org/) to remove duplicates. Due to the structure of the integrated
alignment file the *samtools fixmate -p* option must be disabled. Removing PCR duplicates using paired sequencing
reads will give better results; however, duplicate removal using single end reads can be overly aggressive when
used on bisulfite sequencing data and should be performed on a case by case basis.

```shell
# fixmates to prepare for duplicate removal, use -p to disable proper pair check
samtools fixmate -p -m BSB_pe_test.bam BSB_pe_test.fixmates.bam
# sort bam by coordinates for duplicate calling
samtools sort -@ 2 -o BSB_pe_test.sorted.bam BSB_pe_test.fixmates.bam
# remove duplicate reads
samtools markdup -r BSB_pe_test.sorted.bam BSB_pe_test.dup.bam
# index bam file for methylation calling
samtools index BSB_pe_test.dup.bam
```

### BSBolt CallMethylation

Methylation calling outputs a .CGmap file by default, and makes calls using
all mapped reads not flagged as duplicate at a position by default. 

#### **BSB CallMethylation Commands**

```shell
bsbolt CallMethylation -I {input.bam} -DB {bsbolt DB} -O {output prefix}

-h, --help     show this help message and exit

Input / Output Options:
  -I File         input BAM, input file must be in BAM format with index file
  -DB File        path to index directory
  -O File         output prefix
  -text           output plain text files [False]
  -BG             output calls in bedGraph format [False]
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
```

#### **Methylation Calling**

Methylation calling is performed by counting the number of bisulfite converted bases relative to the number of reads
observed at each cytonsine. Relative to the reference genome methylation status at a cytosine and guanine
can only be called using reads mapped to Watson and Crick strands respectively.

```shell
# Methylation Calling with 2 threads, 
python3 -m bsbolt CallMethylation -I ~/tests/BSB_pe_test.sorted.bam -O ~/tests/BSB_pe_test -DB ~/tests/TestData/BSB_Test_DB -t 2 -verbose > methylation_stats.txt
```

### **Output File**

**CGmap**

CGmap is a tab separated text format describing the methylation status of observed cyotsines

  1. Chromosome
  2. Nucleotide, C for reads mapped to the Watson (sense) strand and G for reads mapped to the Crick (anti-sense) strand
  3. Position, base-pairs from start
  4. Context, three base pair methylation context
  5. Sub-Context, two base pair methylation context
  6. Methylation Value, proportion of methylated bases to total bases
  7. Methylated Bases, methylated nucleotides observed
  8. All Bases, total number of nucleotides observed at the mapping position

```text
chrom   nucleotide  position   context sub-context  methylation_value methylated_bases  all_bases
chr11	G	422436	CHH	CC	0.1	1	10
chr12	G	389290	CHH	CT	0.0	0	10
chr13	G	200552	CHH	CT	0.0	0	10
chr11	C	142826	CG	CG	0.0	0	10
```

**Bedgraph**

BedGraph is functionally similar to CGmap format without cytosine context information and zero indexed positions.

  1. Chromosome
  2. Start Position
  4. End Position
  6. Methylation Percentage, percentage of methylated bases to total observed bases
  7. Methylated Bases, methylated nucleotides observed
  8. Unmethylated Bases, total unmethylated bases

