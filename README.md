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

BSBolt alignment is performed using a modified version of the [BWA-mem2](https://github.com/bwa-mem2/bwa-mem2). Simulation 
relies utilizes [wgsim](https://github.com/lh3/wgsim) to simulate Illumina reads before bisulfite conversion. 


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

*External Dependencies*
- BWA-MEM2 and WGSIM forks compiled using gcc 8.3.1 with zlib 
- Setup.py will fail without gcc and zlib accessible

Python Dependencies
- [pysam](https://github.com/pysam-developers/pysam) >= 0.15.4
- [numpy](https://numpy.org/) >=1.16.3
- [tqdm](https://github.com/tqdm/tqdm) >= 4.31.1

