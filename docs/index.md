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

## Release Notes
- BSBolt v1.4.5
  - Fixed maximum read depth bug that prevented methylation call on site covered by greater than 8000 reads
  - Refactored build script, with experimental support for M1 Macs
- BSBolt v1.4.4
  - The default entry point for BSBolt has changed from **BSBolt** to **bsbolt** for conda compatibility

## **Installation**

### **PyPi Installation**

Pre-compiled binaries can be installed using PyPi. Binaries are provided for python >=3.6
on unix like systems (macOS >=10.15 and linux).

```shell
pip3 install bsbolt --user
```

### **Conda Installation**

BSBolt can be installed using the conda package manager using the installation instructions below. 

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

Installation from source requires xcode command line utilities, [homebrew](https://brew.sh/) macOS package manager, autoconf, and
automake are installed. Xcode through the mac App Store, running the xcode installation command listed below, 
or as part of the [homebrew](https://brew.sh/) macOS package manager installation. The full installation process
can be completed as outlined below.

```shell
# install xcode utilities
xcode-select --install
# install homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# install autoconf
brew install autoconf
# install automake
brew installa automake
# optionally install python
brew install python3.8
# install bsbolt
pip3 install bsbolt
```
## Usage

Following installation BSBolt can be called by invoking the BSBolt module.

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
