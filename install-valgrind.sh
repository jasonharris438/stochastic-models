#! /usr/bin/env bash
# Handles installing the valgrind application to detect memory leaks in C++ programs.

set -euo pipefail
IFS=$'\n\t'

valgrind_version="3.27.1"

wget "https://sourceware.org/pub/valgrind/valgrind-${valgrind_version}.tar.bz2"
valgrind_sha512="4522e345fe31bc10478f21ab261cf7b3e509d80fe98fd0c7c9778ac9a9f90e1d8fa6e54a37c2e23114c01c95131035103fc1e78661710d7f45f5d564a31c1015"

if ! echo "$valgrind_sha512 valgrind-${valgrind_version}.tar.bz2" | sha512sum -c -; then
	echo "Valgrind checksum validation failed."
	exit 1
fi

tar -xjvf "valgrind-${valgrind_version}.tar.bz2"
cd "valgrind-${valgrind_version}" || exit 1
./configure
make
make install
