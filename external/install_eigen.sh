#!/bin/bash

# The tarball and the destination directory
TARFILE="eigen-3.4.0.tar.gz"
DESTDIR="$HOME/opt/"

# Check if the destination directory exists, if not create it
if [ ! -d "$DESTDIR" ]; then
    mkdir -p "$DESTDIR"
fi

# Untar the tarball
tar -xvzf "$TARFILE"

# Install
cd eigen-3.4.0
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=${DESTDIR}/include
make install
