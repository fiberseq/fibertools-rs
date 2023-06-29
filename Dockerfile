FROM ubuntu:20.04

# Update default packages
RUN apt-get update

# stop time prompt 
ARG DEBIAN_FRONTEND=noninteractive

# Get Ubuntu packages
RUN apt-get install -y \
    build-essential \
    curl \
    wget \
    zip \
    cmake \
    git \
    autoconf \
    automake \
    make \
    gcc \
    perl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev

# Update new packages
RUN apt-get update

# Install HTSlib
RUN wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2
RUN tar -xjvf htslib-1.17.tar.bz2
WORKDIR /htslib-1.17
RUN ./configure
RUN make
RUN make install

# Install Rust
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"


# Install libtorch
WORKDIR /src/fibertools-rs
RUN wget https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-2.0.1%2Bcpu.zip
RUN unzip libtorch-shared-with-deps-2.0.1+cpu.zip
ENV LIBTORCH=/src/fibertools-rs/libtorch
ENV LD_LIBRARY_PATH=/src/fibertools-rs/libtorch/lib:$LD_LIBRARY_PATH
ENV LIBTORCH_CXX11_ABI=0

COPY . .

# RUN cargo install --path . --root /usr/local
RUN cargo build --release  
RUN cp ./target/release/ft /usr/local/bin/ft
RUN ft --help