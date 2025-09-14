#### Example of Multistage Build for Various C++ Projects ####
#### Development Environment = builder.
#### 3rd Party Application = installer,
#### Production Environment = Final.

# Parent image ubuntu latest.
FROM ubuntu:latest AS builder

# Install initial dependencies.
COPY setup-env.sh ./
RUN chmod +x setup-env.sh
RUN ./setup-env.sh

# Copy cmake setup script.
COPY install-cmake.sh ./

# Run cmake setup commands.
RUN chmod +x install-cmake.sh
RUN ./install-cmake.sh

# Set workdir.
WORKDIR /usr/src/app

# Copy source code.
COPY src ./src
COPY include ./include
COPY tests ./tests
COPY CMakeLists.txt ./
COPY Config.cmake.in ./

# Build C++ code.
RUN cmake .
RUN cmake --build . --target install

# Final image.
FROM ubuntu:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    libc6 \
    libgsl-dev

# Set workdir.
WORKDIR /usr/src/app

# Copy source code.
COPY --from=builder /usr/src/app/install/include/boost /usr/local/lib
COPY --from=builder /usr/src/app/install/include/stochastic_models /usr/local/lib
COPY --from=builder /usr/src/app/src/libstochastic_models.so /usr/local/lib
