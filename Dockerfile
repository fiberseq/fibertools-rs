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
    git

# Update new packages
RUN apt-get update

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

#RUN cargo install --path .
#FROM mambaorg/micromamba:1.4.4
# mamba setup
#COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
#RUN micromamba install -y -n base -f /tmp/env.yaml && \
#micromamba clean --all --yes 

# build directory setup
#RUN mkdir -p /tmp/fibertools-rs
#RUN chown $MAMBA_USER /tmp/fibertools-rs
#USER $MAMBA_USER
#WORKDIR /tmp/fibertools-rs
#COPY --chown=$MAMBA_USER:$MAMBA_USER  . .

# activate mamba in the build env  
#ARG MAMBA_DOCKERFILE_ACTIVATE=1 

# tell tch-rs to use python libtorch
#ENV LIBTORCH_USE_PYTORCH=1

# test pytorch
#RUN python -c "import torch"

#CMD ["myapp"]