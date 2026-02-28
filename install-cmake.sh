#! /usr/bin/env bash

# Install system dependencies.
# GNU compiler tools, Cmake, the GDB debugger,
# and GNU scientific lib.

set -euo pipefail
IFS=$'\n\t'

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
version="4.2"
build="3"
checksum="7efaccde8c5a6b2968bad6ce0fe60e19b6e10701a12fce948c2bf79bac8a11e9"

# Install cmake.
mkdir ~/temp
cd ~/temp || exit 1

wget "${src_url}/v${version}.${build}/cmake-${version}.${build}.tar.gz"

if ! sha256sum -c <(echo "$checksum cmake-${version}.${build}.tar.gz"); then
	echo "CMake checksum validation failed."
	exit 1
fi

tar -xzvf "cmake-${version}.${build}.tar.gz"
cd "cmake-${version}.${build}/" || exit 1

./bootstrap
make -"j$(nproc)"
make install
