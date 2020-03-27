# BSBolt (BiSulfite Bolt) - Under Construction
## A fast and safe read alignment platform for bisulfite sequencing data

[BiSuflite Bolt (BSBolt)](https://github.com/NuttyLogic/BSBolt); a fast and scalable bisulfite sequencing analysis platform. BSBolt is an integrated analysis 
platform that offers support for bisulfite sequencing read simulation, alignment, methylation calling, data aggregation, 
and data imputation. BSBolt has been validated to work with a wide array of bisulfite sequencing data, 
including whole genome bisulfite sequencing (WGBS), reduced representative bisulfite sequencing data (RRBS), 
and targeted methylation sequencing data. BSBolt utilizes forked versions of [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) 
and [WGSIM](https://github.com/lh3/wgsim) for read alignment and read simulation respectively. BSBolt is release under the 
 MIT licensed. 


### Installation

BSBolt can be installed through the python packages index using the command below.

```shell
pip3 install BSBolt -v
```

BSBolt can also be installed by cloning the github repository and manually compiling included dependencies. 

```shell
# clone the repository
git clone https://github.com/NuttyLogic/BSBolt.git
cd BSBolt
# compile executables and install package
pip3 install .
```

1. Python Dependencies
    - [pysam](https://github.com/pysam-developers/pysam) >= 0.15.2
    - [numpy](https://numpy.org/) >=1.16.3
    - [tqdm](https://github.com/tqdm/tqdm) >= 4.31.1
    
### Usage
Following installation BSBolt can be called by invoking the BSBolt module, 
please see [tutorial][https://bsbolt.readthedocs.io/en/latest/tutorial/].

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
Sort                Sort BAM File
BamIndex            Index BAM file
```




