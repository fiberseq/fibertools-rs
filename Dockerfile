FROM mambaorg/micromamba:1.4.4

# mamba setup
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes 

# build directory setup
RUN mkdir -p /tmp/fibertools-rs
RUN chown $MAMBA_USER /tmp/fibertools-rs
USER $MAMBA_USER
WORKDIR /tmp/fibertools-rs
COPY --chown=$MAMBA_USER:$MAMBA_USER  . .

# activate mamba in the build env  
ARG MAMBA_DOCKERFILE_ACTIVATE=1 

# tell tch-rs to use python libtorch
ARG LIBTORCH_USE_PYTORCH=1

# test pytorch
RUN python -c "import torch"

#RUN cargo build
#RUN cargo clippy --release
#RUN cargo test --release
#CMD ["myapp"]