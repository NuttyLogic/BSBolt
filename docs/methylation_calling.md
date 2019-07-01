## Methylation Calling 

Methylation values, the proportion of methylated bases to the total number of observed bases, are called using the BSBolt CallMethylation module. 
Methylation values for guanine nucleotides are made for reads aligned to the Crick (anti-sense) strands and calls for cytosine nucleotides are made for 
reads aligned to the Watson (sense) strand. Methylation calls for WGBS and targeted bisuflite sequencing can be improved by removing PCR duplicate reads 
before calling methylation. It is not recommend to remove duplicated for RRBS as the sequencing reads will often share the same mapping coordinates 
due to enzymatic digestion. 

### Methylation Calling Pre-Processing

We recommend using [samtools](http://www.htslib.org/) to remove duplicates. Due to the structure of the integrated alignment file the *samtools fixmate -p* option must be disabled. 
Removing PCR duplicates using paired sequencing reads will give better results; removal using single end reads can be overly aggressive when used with bisulfite sequencing data
 and should be performed on a case by case basis. 

```shell
# fixmates to prepare for duplicate removal, use -p to disable proper pair check
samtools samtools fixmate -p -m BSB_pe_test.bam BSB_pe_test.fixmates.bam 
# sort bam by coordinates for duplicate calling
samtools sort -@ 2 -o BSB_pe_test.sorted.bam BSB_pe_test.fixmates.bam
# remove duplicate reads
samtools markdup -r BSB_pe_test.sorted.bam BSB_pe_test.dup.bam
# index bam file for methylation calling
samtools index BSB_pe_test.dup.bam
```

### BSBolt CallMethylation

Methylation calling outputs a .CGmap file by default. To maintain compatibility with some downstream analysis tools 
ATCGmap files can by output, but this feature will be removed in a future update

**BSB CallMethylation Commands**
```shell
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
```shell
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