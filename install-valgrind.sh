#! /usr/bin/env bash
# Handles installing the valgrind application to detect memory leaks in C++ programs.

set -euo pipefail
IFS=$'\n\t'

valgrind_version="3.26.0"

wget "https://sourceware.org/pub/valgrind/valgrind-${valgrind_version}.tar.bz2"
valgrind_md5="856da1bc568212df6df502295a0439c0"

if ! echo "$valgrind_md5 valgrind-${valgrind_version}.tar.bz2" | md5sum -c -; then
	echo "Valgrind checksum validation failed."
	exit 1
fi

tar -xjvf "valgrind-${valgrind_version}.tar.bz2"
cd "valgrind-${valgrind_version}" || exit 1
./configure
make
make install
