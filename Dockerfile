# Parent image gcc:12.
FROM gcc:12

# Install GNU compiler tools, Cmake, the GDB debugger, and GNU scientific lib.
RUN apt-get -yqq update && \
    apt-get -yqq install build-essential gdb libtool autoconf unzip wget &&\
    apt-get -yqq install libgsl-dev

RUN apt remove --purge --auto-remove cmake

# Copy source code.
COPY stochastic_models ./stochastic_models
COPY CMakeLists.txt ./
COPY main.cpp ./
COPY setup.sh ./

# Run cmake setup commands.
RUN chmod +x setup.sh
RUN ./setup.sh

# Build C++ code.
RUN cmake .
RUN cmake --build .

# # Final image.
# FROM ubuntu:23.04 as Final

# # Set workdir.
# WORKDIR /usr/src/app

# # Copy final build.
# COPY --from=builder stochastic-models ./

# CMD [ "./stochastic-models"]
