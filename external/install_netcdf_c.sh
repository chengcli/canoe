#!/bin/bash

# The tarball and the destination directory
TARFILE="netcdf-c-4.9.3.tar.gz"
DESTDIR="$HOME/opt/"

# Check if the destination directory exists, if not create it
if [ ! -d "$DESTDIR" ]; then
    mkdir -p "$DESTDIR"
fi

# Untar the tarball
tar -xvzf "$TARFILE"

# Change to the directory (assuming it's named "netcdf-c-4.9.3" after untar)
cd netcdf-c-4.9.3

# Configure, make and install
FFLAGS="-fPIC -g -O2" FCFLAGS="-fPIC -g -O2" CFLAGS="-fPIC -g -O3" CXXFLAGS="-fPIC -g -O2" ./configure --prefix="$DESTDIR" --disable-hdf5
make
make install

# Print out a message when done
if [ $? -eq 0 ]; then
    echo "netcdf has been successfully installed to $DESTDIR"
else
    echo "There was a problem during the installation."
fi
