FROM gitpod/workspace-full
FROM rust:latest
RUN sudo apt-get install --yes libgsl0-dev

USER gitpod
ENV SHELL /bin/bash
RUN wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
RUN ./bin/micromamba create -n rbt bcftools starcode -c conda-forge -c bioconda
RUN yes | ./bin/micromamba shell init -s bash -p /home/gitpod/micromamba
RUN mamba conda install pytorch==2.0.1 torchvision torchaudio cpuonly -c pytorch 
RUN export LIBTORCH_USE_PYTORCH=1