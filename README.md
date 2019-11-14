# BSBolt (BiSulfite Bolt)
## A fast and safe read alignment platform for bisulfite sequencing data

[BiSuflite Bolt (BSBolt)](https://github.com/NuttyLogic/BSBolt); a fast and scalable bisulfite sequencing analysis platform. BSBolt is an integrated analysis 
platform that offers support for bisulfite sequencing read simulation, alignment, methylation calling, data aggregation, 
and data imputation. BSBolt has been validated to work with a wide array of bisulfite sequencing data, 
including whole genome bisulfite sequencing (WGBS), reduced representative bisulfite sequencing data (RRBS), 
and targeted methylation sequencing data.
 
## Documentation

BSBolt documentation can be found at [bsbolt.readthedocs.io](https://bsbolt.readthedocs.io).

### Installation

At installation two external dependencies 
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
and [ART Read Simulation Tools](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) 
are compiled under the terms of the [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.html).

Install through PIP.

```shell
pip3 install BSBolt -v
```

Install from github.
```shell
git clone git@github.com:NuttyLogic/BSBolt.git
cd BSBolt-master
python3 setup.py install .
```


1. Python Dependencies
    - [pysam](https://github.com/pysam-developers/pysam) >= 0.15.2
    - [numpy](https://numpy.org/) >=1.16.3
    - [tqdm](https://github.com/tqdm/tqdm) >= 4.31.1
2. C Dependencies
    - libgsl
