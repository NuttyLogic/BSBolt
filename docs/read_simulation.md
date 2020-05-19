#### **Simulated Bisulfite Read Meta-data**

Contained in each fastq comment is colon separated meta-data (contig, read start position, read end position,
methyl "cigar", and the reference strand / conversion pattern). The methyl "cigar" contains base level information about
modifiable bases for read reference strand (Watson or Crick).

#### **Methyl Cigar Format**

The methyl cigar is base matched, so unlike a traditional cigar string deletions are not represented. Insertions are
represented numerically, with an integer indicating an inserted base and the position of the base within the insertion.
Sequencing error is simulated after setting methylation, so the methyl cigar is representative of the sequence before
error simulation.

| Operation | Description    | Consumes Query | Consumes Reference |
| :---: | :---:    | :---: | :---: |
|M    |base match|  yes| yes|  
|X    |sequence mismatch|   yes| yes|
|1-9  |insertion and position| yes | no|
|E    |sequence error|  yes| yes|
|C    |methylated CG|   yes| yes|
|c    |unmethylated CG| yes| yes|
|Y    |methylated CH|   yes| yes|
|y    |unmethylated CH| yes| yes|
|Z    |methylated mismatch| yes| yes|
|z    |unmethylated mismatch|   yes| yes|
|R    |methylated insertion|    yes| no|
|r    |unmethylated insertion|  yes| no|
|V    |sequence mismatch that generates false methylation signal (C to T / G to A)|   yes|   yes|


#### **BSB Simulate Commands**

```shell
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
  -NF Float   cutoff threshold for amibiguous bases, simulated reads with a proportion of ambiguous bases above this threshold will not be output [0.05]
  -FM Int     max fragment size [400]
  -IM Int     insert length mean [50]
  -SM Int     insert length standard deviation [50]
```

#### **Simulate Paired End, Undirectional Bisulfite Reads**

```shell
python3 -m BSBolt Simulate -G ~/Tests/TestData/BSB_test.fa -O ~/Tests/TestSimulations/BSB_pe -U -PE
```
