BSBolt consists of mixture of C, CPP, and Python modules. BSBolt was built to be
[manylinux2010](https://www.python.org/dev/peps/pep-0571/) compliant.

### **PyPi Installation**

The easiest installation method is installing pre-compiled binaries using PyPi. Binaries are provided for python >=3.6
on unix like systems (macOS and linux).

```shell
pip3 install BSBolt --user
```

### **Installing from Source**

Dependencies

* zlib-devel >= 1.2.3-29
* GCC >= 8.3.1

```shell
# clone the repository
git clone https://github.com/NuttyLogic/BSBolt.git
cd BSBolt
# compile and install package
pip3 install .
```
