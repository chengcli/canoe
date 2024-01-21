# Canoe: Comprehensive Atmosphere N' Ocean Engine

[![build](https://github.com/chengcli/canoe/actions/workflows/main.yml/badge.svg)](https://github.com/chengcli/canoe/actions/workflows/main.yml)
[![build](https://github.com/chengcli/canoe/actions/workflows/mac.yml/badge.svg)](https://github.com/chengcli/canoe/actions/workflows/mac.yml)
[![DocStatus](https://readthedocs.org/projects/pycanoe/badge/?version=latest)](https://pycanoe.readthedocs.io/en/latest/?badge=latest)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![codecov](https://codecov.io/gh/chengcli/canoe/branch/main/graph/badge.svg?token=hKnnv79a09)](https://codecov.io/gh/chengcli/canoe)

## Install system libraries and toolchain
Canoe can be installed on either a Linux distribution or on MacOS. Open a Linux or Mac terminal,
you can clone this repo using the following command:
```
git clone https://github.com/chengcli/canoe
```
This will copy all source files into your local computer. You will need to install a few
system libraries before installing canoe. All following instructions are executed under
the `canoe/` directory, which is referred to as the `root`.

### MacOS Installation Guide
We assume that [homebrew](https://brew.sh/) is already installed on your Mac and we will
use `brew` to install required system libraries. The system libraries are listed in
`Brewfile` at the `root`. To install them all, execute
```
brew bundle
```
### Ubuntu Linux Installation Guide
On a Ubuntu linux system, use `apt` to install
```
sudo apt install $(cat packages_debian.txt)
```

### Redhat Linux Installation Guide
On a Redhat linux system, use `yum` to install
```
sudo yum -y install $(cat packages_centos.txt)
```

### Multi-core execution
If multi-core parallelization is needed, these extra pacakges should be install
- mpich parallel library

Ubuntu linux:
```
sudo apt install libmpich-dev
```
Redhat linux:
```
sudo yum install mpich-devel
source ~/.bash_profile
```
- pnetcdf output

Redhat linux does not support pnetcdf natively. So it should be downloaded and install.
The default installation directory is $HOME/opt/
```
cd external
./fetch_pnetcdf.sh
./install_pnetcdf.sh
cd ..
```

## Install python libraries
The minimum python version is 3.8.
All needed python libraries are collected in `requirements.txt`. We suggest using a
python [virtual environment](https://docs.python.org/3/library/venv.html) to install
these packages. If you are already using a virtual enviroment, install python packages
by
```
pip3 install -r requirements.txt
```
Otherwise, to create a python virtual environment:
```
python -m venv pyenv
```
This command will create an environment named `pyenv` in your current directory. Then, you
can use the previous command to install the python packages.

## Install pre-commit
Register your `pre-commit` hooks using
```
pre-commit install
```
The contributor's guide explains the meaning of `pre-commit`.

## How to build and test
After you completed the installation steps, you can build the canoe library.
The easiest way is to build it in-place, meaning that the build (binary files) are
located under `root`. To do so, make a new directory named `build`
```
mkdir build
```
All build files will be generated and placed under this directory. It is completely safe
to delete the whole directory if you want another build. `cd` to build and `cmake`

```
cd build
cmake ..
```
This command tells the cmake command to look for `CMakeFiles.txt` in the parent directory,
and start configuring the compile environment. Then compile the code by
```
make -j4
```
This comman will use 4 cores to compile the code in parallel. Once complete, all executable
files will be placed in `build/bin`.

## Optional packages
- The [Reference Forward Model](http://eodg.atm.ox.ac.uk/RFM/) (RFM) is provided optionally as
a tool to generate opacity tables. The source code of this package is not publically available.
Please contact [Anu Dudhia](mailto:anu.dudhia@physics.ox.ac.uk) ar [Cheng Li](mailto:chengcli@umich.edu) to obtain access.
The build process turns off RFM by default. To turn on building RFM, use
```
cmake .. -DRFM=ON
```
- The [DIScrete Ordinate Radiative Transfer](https://doi.org/10.1016/j.jqsrt.2011.03.019) (DISORT)
is provided optionally as a plane-parallel radiative transfer solver.
The original source code was in Fortran77. Tim Downling translated it to C in 2011.
The C-source code, version 2.1.3, is hosted at [libradtran.org](http://libradtran.org/doku.php).
The build process turns off DISORT by default. To turn on building DISORT, use
```
cmake .. -DDISORT=ON
```
- The [Parallel Kernel-Independent Fast Multipole Method](https://github.com/dmalhotra/pvfmm) (PVFMM)
is provided optionally as a potential solver for N-Body simulation.
The upstream source code is hosted at [pvfmm](https://github.com/dmalhotra/pvfmm).
We made a customized fork hosted at [chengcli:pvfmm](https://github.com/chengcli/pvfmm/).
The build process turns off PVFMM by default. To turn on building PVFMM, use
```
cmake .. -DPVFMM=ON
```
