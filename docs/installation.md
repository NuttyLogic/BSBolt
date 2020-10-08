BSBolt and dependencies need to be compiled during the installation process. Pre-compiled binaries are available for 
macOS >= 10.15 and linux distributions released after 2010 according to the 
[manylinux2010](https://www.python.org/dev/peps/pep-0571/) python enhancement proposal. If a precompiled binary is 
unavailable for the target OS binaries will be built from source. If working on a linux distribution zlib-devel should 
be installed prior to compilation. Compilation on macOS requires xcode-command line utilities and autoconf be installed 
as described below.  

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

### **Installing from Source on macOS**

Dependencies 
* autoconf
* homebrew
* xcode

Installation from source requires xcode command line utilities, [homebrew](https://brew.sh/) macOS package manager, 
and autoconf are installed. Xcode through the mac App Store, running the xcode installation command listed below, 
or as part of the [homebrew](https://brew.sh/) macOS package manager installation. The full installation process
can be completed as outlined below.

```shell script
# install xcode utilities
xcode-select --install
# install homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
# install autoconf
brew install autoconf
# optionally install python
brew install python3.8
# install BSBolt
pip3 install BSBolt
```
