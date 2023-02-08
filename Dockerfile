# Parent image gcc:12.
FROM gcc:12

# Install GNU compiler tools, Cmake, and the GDB debugger.
RUN apt-get -yqq update && \
    apt-get -yqq install build-essential gdb libtool autoconf unzip wget

RUN apt remove --purge --auto-remove cmake

# Set workdir.
WORKDIR /usr/src/app

COPY ./setup.sh ./
RUN chmod +x setup.sh

RUN ./setup.sh

CMD ["/bin/bash"]
