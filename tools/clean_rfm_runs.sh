#! /bin/bash

find . -mindepth 1 -maxdepth 1 -type d -exec rm -r {} \;
find . -name \*.hit -exec rm {} \;
find . -name rfm.log -exec rm {} \;
find . -name rfm.drv -exec rm {} \;
find . -name rfm.atm -exec rm {} \;
find . -name \*.in -exec rm {} \;
