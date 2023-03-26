# Parent image gcc:12.
FROM gcc:12

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
