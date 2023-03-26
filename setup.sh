#! /usr/bin/env bash

# Install system dependencies.
# GNU compiler tools, Cmake, the GDB debugger,
# and GNU scientific lib.
apt-get -yqq update && \
    apt-get -yqq install build-essential gdb libtool autoconf unzip wget &&\
    apt-get -yqq install libgsl-dev

# Remove any existing cmake installation.
apt remove --purge --auto-remove cmake

# Set cmake configuration variables.
version=3.26
build=1

# Install cmake.
mkdir ~/temp
cd ~/temp
wget https://cmake.org/files/v$version/cmake-$version.$build.tar.gz
tar -xzvf cmake-$version.$build.tar.gz
cd cmake-$version.$build/

./bootstrap
make -j$(nproc)
make install
