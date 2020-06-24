# **BSBolt (BiSulfite Bolt)**
## A fast and safe bisulfite sequencing processing platform

[BiSuflite Bolt (BSBolt)](https://github.com/NuttyLogic/BSBolt); a fast and scalable bisulfite sequencing 
processing platform. BSBolt is an integrated analysis platform that offers support for bisulfite sequencing 
read simulation, alignment, methylation calling, data aggregation, and data imputation. BSBolt has been validated 
to work with a wide array of bisulfite sequencing data,including whole genome bisulfite sequencing (WGBS), 
reduced representative bisulfite sequencing data (RRBS), and targeted methylation sequencing data. 
BSBolt utilizes forked versions of [BWA](https://github.com/lh3/bwa) 
and [WGSIM](https://github.com/lh3/wgsim) for read alignment and read simulation respectively. 
BSBolt is released under the MIT license.
 
## Documentation

BSBolt documentation can be found at [bsbolt.readthedocs.io](https://bsbolt.readthedocs.io).

### Installation

BSBolt alignment is performed using a modified version of [BWA](https://github.com/lh3/bwa). Read simulation is carried 
out using a modified version of [wgsim](https://github.com/lh3/wgsim). 


Install through PyPi

```shell
pip3 install BSBolt --user
```

Install from source
```shell
git clone git@github.com:NuttyLogic/BSBolt.git
cd BSBolt-master
python3 setup.py install .
```

