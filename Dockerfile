FROM ubuntu:22.04

ENV FT_VERSION="0.3.6"

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
# RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2
# RUN tar -xjvf htslib-1.18.tar.bz2
# WORKDIR /htslib-1.18
# RUN ./configure
# RUN make
# RUN make install

# Install Rust
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"


# Install libtorch
WORKDIR /src/fibertools-rs
ENV TORCH_VERSION="2.1.0"
# "cu117" or "cu118" or "cpu"
ENV INSTALL_TYPE="cu118" 
ENV INSTALL_TYPE="cpu" 
ENV PY_URL="https://download.pytorch.org/libtorch/${INSTALL_TYPE}/libtorch-shared-with-deps-${TORCH_VERSION}%2B${INSTALL_TYPE}.zip"
RUN wget ${PY_URL}
RUN unzip "libtorch-shared-with-deps-${TORCH_VERSION}+${INSTALL_TYPE}.zip"
RUN ls
#RUN unzip libtorch*.zip
ENV LIBTORCH=/src/fibertools-rs/libtorch
ENV LD_LIBRARY_PATH=/src/fibertools-rs/libtorch/lib:$LD_LIBRARY_PATH
ENV LIBTORCH_CXX11_ABI=0

COPY . .

#RUN wget "https://github.com/fiberseq/fibertools-rs/archive/refs/tags/v${FT_VERSION}.tar.gz"
#RUN tar -xzf "v${FT_VERSION}.tar.gz"
#RUN cd fibertools-rs-${FT_VERSION} && ls && 
RUN cargo build --release && cp ./target/release/ft /usr/local/bin/ft

RUN ft --help

#RUN cargo install --path . --root /usr/local