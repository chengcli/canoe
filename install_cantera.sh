#! /usr/bin/bash

CXX=g++
CC=gcc
cxx_flags=-std=c++14
prefix=${HOME}/opt/
python_package=n
f90_interface=n
system_eigen=y
extra_inc_dirs=${HOME}/opt/include/eigen3
system_blas_lapack=n
boost_inc_dir=`pwd`/cliboost/

cd cantera
scons build CXX=${CXX} CC=${CC} cxx_flags=${cxx_flags} prefix=${prefix} \
  python_package=${python_package} f90_interface=${f90_interface} \
  system_eigen=${system_eigen} extra_inc_dirs=${extra_inc_dirs} \
  system_blas_lapack=${system_blas_lapack} \
  boost_inc_dir=${boost_inc_dir}
scons install
cd ..
