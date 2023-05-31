#! /bin/bash

find . -mindepth 1 -maxdepth 1 -type d -exec rm -r {} \;
find . -name \*.hit -exec rm {} \;
find . -name \*.log -exec rm {} \;
find . -name \*.inp -exec rm {} \;
find . -name \*.inp -exec rm {} \;
find . -name kcoeff.inp-\* -exec rm {} \;
