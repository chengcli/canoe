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

# harp
cpplint --filter=${filters} --recursive src/inversion
