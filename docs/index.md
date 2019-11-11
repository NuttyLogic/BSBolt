# BSBolt (BiSulfite Bolt) - Under Construction
## A fast and safe read alignment platform for bisulfite sequencing data

[BiSuflite Bolt (BSBolt)](https://github.com/NuttyLogic/BSBolt); a fast and scalable bisulfite sequencing analysis platform. BSBolt is an integrated analysis 
platform that offers support for bisulfite sequencing read simulation, alignment, methylation calling, data aggregation, 
and data imputation. BSBolt has been validated to work with a wide array of bisulfite sequencing data, 
including whole genome bisulfite sequencing (WGBS), reduced representative bisulfite sequencing data (RRBS), 
and targeted methylation sequencing data. BSBolt utilizes 
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) as its alignment engine and 
[ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) for simulation of 
Illumina reads. 

### Installation

At installation BSBolt compiles two external dependencies 
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
and [ART Read Simulation Tools](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) 
under the terms of the [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.html). 

```shell
pip3 install BSBolt -v
```

1. Python Dependencies
    - [pysam](https://github.com/pysam-developers/pysam) >= 0.15.2
    - [numpy](https://numpy.org/) >=1.16.3
    - [tqdm](https://github.com/tqdm/tqdm) >= 4.31.1
2. External C Dependencies
    - libgsl
    - zlib      
    - libbz2
    - liblzma
    - libcurl
    - libcrypto  
 
### Usage
Following installation BSBolt can be called by invoking the BSBolt module. Please see the [tutorial][https://bsbolt.readthedocs.io/en/latest/tutorial/] for 
a detailed walk through.

```shell
python3 -m BSBolt
```
 
```shell
BSBolt Module

Align               Alignment Module
Index               Index Generation Module
CallMethylation     Methylation Calling Module
AggregateMatrix     CGmap Matrix Aggregation Module
Simulate            BSBolt Illumina Read Simulation Module
Impute              kNN Imputation Module
```




