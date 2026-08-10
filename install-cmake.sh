#! /usr/bin/env bash

# Install system dependencies and the prebuilt CMake binary.
# GNU scientific lib, git for FetchContent, CMake from the
# official release archive.

set -euo pipefail
IFS=$'\n\t'

apt-get -yqq update && apt-get -yqq install \
	git \
	ca-certificates \
	wget

# Must follow the first install block.
apt-get -yqq install libgsl-dev

# Remove any existing cmake installation.
apt-get remove --purge --auto-remove cmake

# Set cmake configuration.
src_url="https://github.com/Kitware/CMake/releases/download"
version="4.2"
build="3"
checksum="5bb505d5e0cca0480a330f7f27ccf52c2b8b5214c5bba97df08899f5ef650c23"
archive="cmake-${version}.${build}-linux-x86_64.tar.gz"

# Install the prebuilt cmake binary.
mkdir ~/temp
cd ~/temp || exit 1

wget "${src_url}/v${version}.${build}/${archive}"

if ! sha256sum -c <(echo "$checksum $archive"); then
	echo "CMake checksum validation failed."
	exit 1
fi

tar -xzf "$archive" --strip-components=1 -C /usr/local
