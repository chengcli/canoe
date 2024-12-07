#! /bin/bash

CXX=g++
CC=gcc
prefix=${HOME}/opt/
boost_inc_dir=`pwd`/cantera/ext/cliboost/
python_package=n
f90_interface=n
system_eigen=n
system_blas_lapack=n
hdf_support=n
python_sdist=n

cd cantera
git submodule init ext/cliboost
git submodule update ext/cliboost

scons build CXX=${CXX} CC=${CC} cxx_flags="-std=c++17 -D_GLIBCXX_USE_CXX11_ABI=0" \
  python_package=${python_package} \
  python_sdist=${python_sdist} \
  f90_interface=${f90_interface} \
  system_eigen=${system_eigen} \
  system_blas_lapack=${system_blas_lapack} \
  hdf_support=${hdf_support} \
  prefix="${prefix}" \
  boost_inc_dir="${boost_inc_dir}" -j8
scons install
cd ..
