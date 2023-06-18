FROM mambaorg/micromamba:1.4.4

# mamba setup
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes 

# activate mamba in the build env  
ARG MAMBA_DOCKERFILE_ACTIVATE=1 
RUN python -c 'import uuid; print(uuid.uuid4())' > /tmp/my_uuid

# env setup 
RUN echo "LIBTORCH"
RUN echo $LIBTORCH_USE_PYTORCH
RUN export LIBTORCH_USE_PYTORCH=1
RUN echo $LIBTORCH_USE_PYTORCH
RUN echo "END LIBTORCH"

# rust setup
WORKDIR /usr/src/myapp
COPY . .

RUN cargo install --path .

CMD ["myapp"]