#### Multistage build: the builder compiles and installs the library,   ####
#### the final image ships the install prefix as a build-against base.  ####

FROM ubuntu:24.04@sha256:561618e2c15bf2397621dd04f96926663a3b5616c189cf7e38db7e82f5c538ea AS builder

# Install the gcc toolchain.
COPY setup-env.sh ./
RUN chmod +x setup-env.sh
RUN ./setup-env.sh

# Install the build tools and GSL.
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    ca-certificates \
    libgsl-dev

# Set workdir.
WORKDIR /usr/src/app

# Copy source code.
COPY src ./src
COPY include ./include
COPY examples ./examples
COPY CMakeLists.txt ./
COPY Config.cmake.in ./

# Build and stage the install tree.
RUN cmake -S . -B build -DBUILD_TESTING=OFF
RUN cmake --build build -j"$(nproc)"
RUN cmake --install build --prefix /opt/stochastic-models

# Final image.
FROM ubuntu:24.04@sha256:561618e2c15bf2397621dd04f96926663a3b5616c189cf7e38db7e82f5c538ea

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

# Ship the install prefix.
COPY --from=builder /opt/stochastic-models/include /usr/local/include
COPY --from=builder /opt/stochastic-models/lib /usr/local/lib

# Register the library chain.
RUN ldconfig

USER ubuntu
WORKDIR /home/ubuntu
