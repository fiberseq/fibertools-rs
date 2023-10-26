FROM pytorch/pytorch:2.1.0-cuda12.1-cudnn8-runtime

# stop time prompt 
ARG DEBIAN_FRONTEND=noninteractive

# Update new packages
RUN apt-get update

# Get Ubuntu packages
RUN apt-get install -y build-essential curl 

# Update new packages
RUN apt-get update

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sed 's#/proc/self/exe#\/bin\/sh#g' | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN cargo --help

# Install fibertools-rs
WORKDIR /src/fibertools-rs
COPY . .

# use the install pytroch from docker image
ENV LIBTORCH_USE_PYTORCH=1
#RUN python checkpy.py
#RUN ls "/opt/conda/lib/python3.10/site-packages/torch/lib/"
#ENV LIBTORCH="/opt/conda/lib/python3.10/site-packages/torch/"
#ENV LD_LIBRARY_PATH="${LIBTORCH}/lib:${LD_LIBRARY_PATH}"

RUN cargo build --release && cp ./target/release/ft /usr/local/bin/ft
RUN ft --help