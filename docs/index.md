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
BSBolt requires [python3.6](python.org) or greater to run, with the external *pysam*, *numpy*, and *tqdm* packages. BSBolt is 
distributed with [Bowtie2]( http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
and [ART Read Simulation Tools](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) 
under the terms of the [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.html). BSBolt is implemented for Unix like
operating systems (MacOS, Linux). Windows support under the Windows Subsystem for Linux (WSL) if functional, but has been 
minimally tested. After installing python >= 3.6 follow the instructions below to complete the installation process.

1. **Clone BSBolt**
    ```shell
    git clone  https://github.com/NuttyLogic/BSBolt.git
    ```
2. **Installing Dependencies**
    ```shell
    cd BSBolt
    pip3 install -r requirements.txt 
    ```

### Usage
BSBolt modules can be invoked  by calling *BSBolt.py*. Module documentation can be called by calling the module followed by 
**-h**.
 
```shell
python3 BSBolt.py Module

Align               Alignment Module
Index               Index Generation Module
CallMethylation     Methylation Calling Module
AggregateMatrix     CGmap Matrix Aggregation Module
Simulate            BSBolt Illumina Read Simulation Module
Impute              kNN Imputation Module
```
