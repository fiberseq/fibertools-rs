FROM pytorch/pytorch:2.1.0-cuda12.1-cudnn8-runtime

# stop time prompt 
ARG DEBIAN_FRONTEND=noninteractive

# Update new packages
RUN apt-get update

# Get Ubuntu packages
RUN apt-get install -y \
    build-essential \
    curl \
    zip \
    cmake \
    git \
    autoconf \
    automake \
    make \
    gcc \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev

# Update new packages
RUN apt-get update

# Install Rust
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

COPY . .

# use the install pytroch from docker image
ENV LIBTORCH_USE_PYTORCH=1
RUN python checkpy.py
RUN ls "/opt/conda/lib/python3.10/site-packages/torch/lib/"

RUN cargo build --release && cp ./target/release/ft /usr/local/bin/ft
RUN ft --help