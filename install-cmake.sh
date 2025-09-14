#! /usr/bin/env bash

# Install system dependencies.
# GNU compiler tools, Cmake, the GDB debugger,
# and GNU scientific lib.
apt-get -yqq update && apt-get -yqq install \
    git \
    gdb \
    libtool \
    autoconf \
    unzip \
    wget

# Must follow the first install block.
apt-get -yqq install libgsl-dev

# Remove any existing cmake installation.
apt remove --purge --auto-remove cmake

# Set cmake configuration.
src_url="https://github.com/Kitware/CMake/releases/download"
version="4.1"
build="1"
checksum="b29f6f19733aa224b7763507a108a427ed48c688e1faf22b29c44e1c30549282"

# Install cmake.
mkdir ~/temp
cd ~/temp || exit 1

wget "${src_url}/v${version}.${build}/cmake-${version}.${build}.tar.gz"

if ! sha256sum -c <(echo "$checksum cmake-${version}.${build}.tar.gz")
then
    echo "CMake checksum validation failed."
    exit 1
fi

tar -xzvf "cmake-${version}.${build}.tar.gz"
cd "cmake-${version}.${build}/" || exit 1

./bootstrap
make -"j$(nproc)"
make install
