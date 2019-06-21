# BSBolt (BiSulfite Bolt) - Under Construction
## A fast and safe read alignment platform for bisulfite sequencing data

[BiSuflite Bolt (BSBolt)](https://github.com/NuttyLogic/BSBolt); a fast and scalable bisulfite sequencing analysis platform. BSBolt is an integrated analysis 
platform that offers support for bisulfite sequencing read simulation, alignment, methylation calling, data aggregation, 
and data imputation. BSBolt has been validated to work with a wide array of bisulfite sequencing data, 
including whole genome bisulfite sequencing (WGBS), reduced representative bisulfite sequencing data (RRBS), 
and targeted methylation sequencing data. 

## Overview
BSBolt provides support for processing and simulation of *Whole Genome Bisulfite Sequencing* (WGBS), *Reduced 
Representative Bisulfite Sequencing* (RRBS), and targeted bisulfite sequencing data. BSBolt utilizes
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) as its alignment engine and 
[ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) for simulation of 
Illumina reads. 

### Installation
BSBolt requires [python3.6](python.org) or greater to run, with the external *pysam*, *numpy*, and *tqdm* packages. BSBolt is 
distributed with [Bowtie2]( http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
and [ART Read Simulation Tools](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) 
under the terms of the [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.html). BSBolt is implemented for Unix like
operating systems (MacOS, Linux). Windows support under the Windows Subsystem for Linux (WSL) works, but has been 
minimally tested. 

**Clone BSBolt**
```shell
git clone  git@github.com:NuttyLogic/BSBolt.git
```

**Installing Dependencies**
```shell
cd BSBolt
pip3 install --requirement requirements.txt 
```

**Calling BSBolt Modules**
Invoke BSBolt modules by calling *BSBolt.py* and selecting the appropriate module. Module help command can be invoked 
by calling *-h* with each module. 
```shell
python3 BSBolt.py Module

Align               Alignment Module
Index               Index Generation Module
CallMethylation     Methylation Calling Module
AggregateMatrix     CGmap Matrix Aggregation Module
Simulate            BSBolt Illumina Read Simulation Module
Impute              kNN Imputation Module
```
