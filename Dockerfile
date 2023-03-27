#### Example of Multistage Build for Various C++ Projects ####
#### Development Environment = builder.
#### 3rd Party Application = installer,
#### Production Environment = Final.

# Parent image gcc:12.
FROM gcc:12 as builder

# Copy source code.
COPY stochastic_models ./stochastic_models
COPY CMakeLists.txt ./
COPY Config.cmake.in ./
COPY setup.sh ./

# Run cmake setup commands.
RUN chmod +x setup.sh
RUN ./setup.sh

# Build C++ code.
RUN cmake .
RUN cmake --build . --target install

# Parent image gcc:12.
FROM gcc:12 as installer

# Copy source code.
COPY --from=builder install ./install
COPY setup.sh ./
COPY subsequent-application/CMakeLists.txt ./
COPY subsequent-application/main.cpp ./

# Run cmake setup commands.
RUN chmod +x setup.sh
RUN ./setup.sh

# Build C++ code.
RUN cmake .
RUN cmake --build .

# Final image.
FROM ubuntu:23.04 as Final

# Set workdir.
WORKDIR /usr/src/app

# Copy final build.
COPY --from=installer another-application ./

CMD [ "./another-application" ]
