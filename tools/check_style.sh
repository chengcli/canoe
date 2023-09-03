#! /bin/bash

filters=-legal/copyright,-build/include_subdir,-build/include_order,-runtime/references

# canoe
cpplint --filter=${filters} src/*.hpp

# astro
cpplint --filter=${filters} --recursive src/astro

# debugger
cpplint --filter=${filters} --recursive src/debugger

# outputs
for file in $(find src/outputs -name '*.cpp' -or -name '*.hpp' | grep -v "mppnccombine.cpp"); do
    cpplint --filter=${filters} $file
done

# utils
for file in $(find src/utils -name '*.cpp' -or -name '*.hpp' | grep -v "ndarrays.hpp"); do
    cpplint --filter=${filters} $file
done

# harp
cpplint --filter=${filters} --recursive src/harp

# opacity
cpplint --filter=${filters} src/opacity/*.cpp

# inversion
cpplint --filter=${filters} --recursive src/inversion

# snap
cpplint --filter=${filters} src/snap/*.hpp
cpplint --filter=${filters} src/snap/*.cpp
cpplint --filter=${filters} --recursive src/snap/eos
cpplint --filter=${filters} --recursive src/snap/thermodynamics
cpplint --filter=${filters} --recursive src/snap/reconstruct
cpplint --filter=${filters} --recursive src/snap/decomposition
cpplint --filter=${filters} --recursive src/snap/implicit
cpplint --filter=${filters} --recursive src/snap/turbulence

# tools
cpplint --filter=${filters} tools/*.cpp

# tracer
cpplint --filter=${filters} --recursive tracer

# c3m
cpplint --filter=${filters} --recursive c3m

# transport
cpplint --filter=${filters} --recursive src/transport

# examples
cpplint --filter=${filters},-runtime/references --recursive examples/1d-rad-jupiter/*.cpp

# tests
cpplint --filter=${filters} tests/test_yaml_read.cpp
