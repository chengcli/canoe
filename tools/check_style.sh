#! /bin/bash

filters=-legal/copyright,-build/include_subdir,-build/include_order
cpplint --filter=${filters} --recursive src/astro
