FROM ubuntu:20.04
FROM rust:latest
FROM mambaorg/micromamba:1.4.4
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

RUN mamba conda install pytorch==2.0.1 torchvision torchaudio cpuonly -c pytorch 
RUN export LIBTORCH_USE_PYTORCH=1

RUN cargo --help