#! /bin/bash

CXX=g++
CC=gcc
cxx_flags="-std=c++17 -D_GLIBCXX_USE_CXX11_ABI=0"
prefix="${HOME}/opt/"
python_package=y
f90_interface=n
system_eigen=n
system_blas_lapack=n
boost_inc_dir=`pwd`/cantera/ext/cliboost/

cd cantera
git submodule init ext/cliboost
git submodule update ext/cliboost

scons build CXX=${CXX} CC=${CC} cxx_flags=${cxx_flags} prefix=${prefix} \
  python_package=${python_package} f90_interface=${f90_interface} \
  system_eigen=${system_eigen} \
  system_blas_lapack=${system_blas_lapack} \
  boost_inc_dir="${boost_inc_dir}" -j8
scons install
cd ..
