#!/bin/bash

# The tarball and the destination directory
TARFILE="pnetcdf-1.12.3.tar.gz"
DESTDIR="$HOME/opt/"

# Check if the destination directory exists, if not create it
if [ ! -d "$DESTDIR" ]; then
    mkdir -p "$DESTDIR"
fi

# Untar the tarball
tar -xvzf "$TARFILE"

# Change to the directory (assuming it's named "pnetcdf-1.12.3" after untar)
cd pnetcdf-1.12.3

NC_INC=${NETCDF_ROOT}/include
NC_LIB=${NETCDF_ROOT}/lib

# Configure, make and install
#FFLAGS="-fPIC -g -O2" FCFLAGS="-fPIC -g -O2" CFLAGS="-fPIC -g -O3" CXXFLAGS="-fPIC -g -O2" ./configure --prefix="$DESTDIR"
CPPFLAGS="-I${MPI_ROOT}/include" LDFLAGS="-L${MPI_ROOT}/lib" ./configure --prefix="$DESTDIR" --with-netcdf4=${NC_ROOT} --with-pic
make -j4
make install

# Print out a message when done
if [ $? -eq 0 ]; then
    echo "pnetcdf has been successfully installed to $DESTDIR"
else
    echo "There was a problem during the installation."
fi
