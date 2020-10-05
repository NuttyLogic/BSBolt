## Bisulfite Sequencing Simulation

Simulation of bisuflite sequencing reads occurs in a series of distinct steps. 

1. Genetic variants are randomly set for each passed contig, unless the mutation rate is set to zero.
    - Heterozygous variants are enable by default, only homozygous variants will be simulated in haplotype mode *-HA*
2. Methylation values are set for modifiable bases (C & G), either randomly or using a provided methylation reference.
    - Random methylation values are set randomly a distinct binomial distribution for CpG and CH sites. *Note*, correlation 
    structure between neighboring CpG sites is not considered when simulating random methylation values.
    - If a reference methylation file is provided as either as a CGmap file or a previously generated BSBolt reference, the 
    methylation value is set by the reference. Simulated methylation of genetic variants is still randomly set.  
3. Reads are randomly generated across the passed contigs by sampling from a uniform distribution to get the read start position. 
4. Sequencing errors are introduced as represented by the base quality in the resulting fastq file.
    - Sequencing errors are not considered modifiable bases during read simulation.
5. Methylation status of modifiable bases is set according to the set methylation value, with the probability of the 
   base being methylation equal to the methylation value.
6. Reads are output as fastq files with individual read meta-data in the fastq comment line.   

**Simulated Bisulfite Read Meta-data**

Contained in each fastq comment is colon separated meta-data (contig, read start position, read end position, 
methyl "cigar", and the reference strand / conversion pattern). The methyl "cigar" contains base level information about 
modifiable bases for read reference strand (Watson or Crick).
        
**Methyl Cigar Format**

The methyl cigar is base matched, so unlike a traditional cigar string deletions are not represented. Insertions are 
represented numerically, with an integer indicating an inserted base and the position of the base within the insertion. 
Sequencing error is simulated after setting methylation, so the methyl cigar is representative of the sequence before 
error simulation.

| Operation | Description    | Consumes Query | Consumes Reference |     
| :---: | :---:    | :---: | :---: | 
|M    |base match|  yes| yes|  
|X    |sequence mismatch|   yes| yes|
|1-9  |Insertion and position| yes | no|    
|E    |sequence error|  yes| yes|
|C    |methylated CG|   yes| yes|
|c    |unmethylated CG| yes| yes|
|Y    |methylated CH|   yes| yes|
|y    |unmethylated CH| yes| yes| 
|Z    |methylated mismatch| yes| yes|
|z    |unmethylated mismatch|   yes| yes|
|R    |methylated insertion|    yes| no|
|r    |unmethylated insertion|  yes| no|

 

**BSB Simulate Commands**
```shell
  -h, --help  show this help message and exit
  -G          Path for reference genome fasta file, fasta file should contain all contigs
  -O          Output prefix
  -PE         Simulate Paired End Reads, default Single End
  -RL         Simulated Read Length
  -RD         Simulated Read Depth
  -U          Simulate Undirectional Reads, default=Directional
  -CG         Path to CGmap file to generate simulation reference profile
  -BR         Path to previously generated BSBolt methylation reference
  -MR         Mutation rate
  -MI         Mutation indel fraction
  -ME         Mutation indel extension probability
  -RS         Random seed for variant generation
  -HA         Haplotype mode, homozygous variants only
  -CH         Skip simulation of CH methylation, all CH sites unmethylated
  -NS         By default observed methylation counts are saved, disable this behavior
  -SE         Sequencing Error
  -NF         Cutoff threshold for amibiguous bases, simulated reads with a proportion of ambiguous bases above this 
              threshold will not be output
  -FM         Max fragment size
  -IM         Insert length mean
  -SM         Insert length standard deviation
  -verbose    Verbose read simulation
  -overwrite  Overwrite previously generated simulation database
```

**Simulate Paired End, Undirectional Bisulfite Reads**
```shell
python3 -m BSBolt Simulate -G ~/tests/TestData/BSB_test.fa -O ~/tests/TestSimulations/BSB_pe -U -PE
```