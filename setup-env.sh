#! /usr/bin/env bash

# Install packages required by gcc.
apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    make \
    ca-certificates \
    libc6 \
    libssl-dev \
    gcc-14 \
    g++-14

# Update default links.
update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100
update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100
