BSBolt and dependencies need to be compiled during the installation process. Pre-compiled binaries are available for 
macOS >= 10.15 and linux distributions released after 2010 according to the 
[manylinux2010](https://www.python.org/dev/peps/pep-0571/) python enhancement proposal. If a precompiled binary is 
unavailable for the target OS binaries will be built from source. Compilation on macOS requires xcode-command line utilities, autoconf, and automake be installed 
as described below.  

### **PyPi Installation**

The easiest installation method is installing pre-compiled binaries using PyPi. Binaries are provided for python >=3.6
on unix like systems (macOS and linux).

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
* autoconf
* automake
* homebrew
* xcode

Installation from source requires xcode command line utilities, [homebrew](https://brew.sh/) macOS package manager, autoconf, and automake are installed. Xcode through the mac App Store, running the xcode installation command listed below, 
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
brew install automake
# optionally install python
brew install python3.8
# install bsbolt
pip3 install bsbolt
```
