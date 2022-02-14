# **BSBolt (BiSulfite Bolt)**

## A fast and safe bisulfite sequencing analysis platform

[BiSuflite Bolt (BSBolt)](https://github.com/NuttyLogic/BSBolt); a fast and scalable bisulfite sequencing
analysis platform. BSBolt is an integrated analysis platform that offers support for bisulfite sequencing
read simulation, alignment, methylation calling, data aggregation, and data imputation. BSBolt has been validated
to work with a wide array of bisulfite sequencing data,including whole genome bisulfite sequencing (WGBS),
reduced representative bisulfite sequencing data (RRBS), and targeted methylation sequencing data.
BSBolt utilizes forked versions of [BWA](https://github.com/lh3/bwa)
and [WGSIM](https://github.com/lh3/wgsim) for read alignment and read simulation respectively.
BSBolt is released under the MIT license.

## Publication

[Farrell, C., Thompson, M., Tosevska, A., Oyetunde, A. & Pellegrini, M. 
**BiSulfite Bolt: A BiSulfite Sequencing Analysis Platform.** 2020.10.06.328559 (2020). 
doi:10.1101/2020.10.06.328559](https://academic.oup.com/gigascience/article/10/5/giab033/6272610)

## Documentation

Documentation can be found at [https://bsbolt.readthedocs.io](https://bsbolt.readthedocs.io/en/latest/).

## Release Notes
- v1.5.0
  - Improved thread handling for methylation / variant calling.
  - Experimental bisulfite aware SNP caller. 
- v1.4.8
  - Fixed bug ending alignment when the reference template end greater than reference boundary.
- v1.4.7
  - Alignment stats fix.
- v1.4.6
  - Alignment statistics now output as generated.
  - Fixed bug where alignment would stop when observed mappability was low.
- v1.4.5
  - Fixed maximum read depth bug that prevented methylation call on site covered by greater than 8000 reads
  - Refactored build script, with experimental support for M1 Macs
- v1.4.4
  - The default entry point for BSBolt has changed from **BSBolt** to **bsbolt** for conda compatibility

## **Installation**

### **PyPi Installation**

Pre-compiled binaries can be installed using PyPi. Binaries are available for python >=3.6
on unix like systems (macOS >=10.15 and linux).

```shell
pip3 install bsbolt --user
```

### **Conda Installation**

BSBolt can be installed using the conda package manager using the instructions below.

```shell
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c cpfarrell bsbolt
```

### **Installing from Source**

Dependencies

* zlib-devel >= 1.2.3-29
* GCC >= 8.3.1

```shell
# clone the repository
git clone https://github.com/NuttyLogic/BSBolt.git
cd bsbolt
# compile and install package
pip3 install .
```

### **Installing from Source on macOS**

Dependencies 

- autoconf
- automake
- homebrew
- xcode

Installation from source requires xcode command line utilities, [homebrew](https://brew.sh/) macOS package manager, 
autoconf, python (>=3.6), and automake.The full installation process is outlined below.

```shell
# install xcode utilities
xcode-select --install
# install homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# install autoconf
brew install autoconf
# install automake
brew install automake
# optionally install python > 3.5
brew install python3.8
# clone the repository
git clone https://github.com/NuttyLogic/BSBolt.git
cd BSBolt
# compile and install package
pip3 install -e .
```
## Usage

Following installation BSBolt can be called using **bsbolt Module**.

```shell
python3 -m bsbolt
```

```shell
bsbolt Module

Align               Alignment
Index               Index Generation
CallMethylation     Methylation Calling
AggregateMatrix     CGmap Matrix Aggregation
Simulate            bsbolt Illumina Read Simulation
Impute              kNN Imputation
Sort                Sort BAM File
BamIndex            Index BAM file
```
