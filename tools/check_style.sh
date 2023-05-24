#! /bin/bash

filters=-legal/copyright,-build/include_subdir,-build/include_order,-runtime/printf
cpplint --filter=${filters} --recursive src/astro
cpplint --filter=${filters} --recursive src/debugger

for file in $(find src/outputs -name '*.cpp' -or -name '*.hpp' | grep -v "mppnccombine.cpp"); do
    cpplint --filter=${filters} $file
done
