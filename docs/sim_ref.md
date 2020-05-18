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