#! /bin/bash

filters=-legal/copyright,-build/include_subdir,-build/include_order,-runtime/printf

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

# inversion
cpplint --filter=${filters} --recursive src/inversion

# snap
cpplint --filter=${filters} src/snap/*.hpp
cpplint --filter=${filters} src/snap/*.cpp
cpplint --filter=${filters} --recursive src/snap/eos
cpplint --filter=${filters} --recursive src/snap/thermodynamics
#cpplint --filter=${filters} --recursive src/snap/reconstruct
#cpplint --filter=${filters} --recursive src/snap/decomposition
#cpplint --filter=${filters} --recursive src/snap/implicit

# tools
cpplint --filter=${filters} tools/*.cpp

# transport
cpplint --filter=${filters} --recursive src/transport

# tests
cpplint --filter=${filters},-runtime/references --recursive tests/1d-rad-jupiter/*.cpp
