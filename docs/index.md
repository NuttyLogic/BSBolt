# **BSBolt (BiSulfite Bolt)**

## A fast and safe bisulfite sequencing processing platform

[BiSuflite Bolt (BSBolt)](https://github.com/NuttyLogic/BSBolt); a fast and scalable bisulfite sequencing
analysis platform. BSBolt is an integrated analysis platform that offers support for bisulfite sequencing
read simulation, alignment, methylation calling, data aggregation, and data imputation. BSBolt has been validated
to work with a wide array of bisulfite sequencing data,including whole genome bisulfite sequencing (WGBS),
reduced representative bisulfite sequencing data (RRBS), and targeted methylation sequencing data.
BSBolt utilizes forked versions of [BWA](https://github.com/lh3/bwa)
and [WGSIM](https://github.com/lh3/wgsim) for read alignment and read simulation respectively.
BSBolt is released under the MIT license.

### Installation

Pre-compiled binaries can be installed through the python package index. Please see the [installation](installation.md)
for instructions about installing from source and package dependency information.

```shell
pip3 install BSBolt -v
```

### Usage

Following installation BSBolt can be called by invoking the BSBolt module.

```shell
python3 -m BSBolt
```

```shell
BSBolt Module

Align               Alignment
Index               Index Generation
CallMethylation     Methylation Calling
AggregateMatrix     CGmap Matrix Aggregation
Simulate            BSBolt Illumina Read Simulation
Impute              kNN Imputation
Sort                Sort BAM File
BamIndex            Index BAM file
```
