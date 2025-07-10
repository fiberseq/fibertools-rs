#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate pytorch2.6
export LIBTORCH_USE_PYTORCH=1
export DYLD_LIBRARY_PATH=/Users/mrvollger/miniconda3/envs/pytorch2.6/lib/python3.13/site-packages/torch/lib/:$DYLD_LIBRARY_PATH
/Users/mrvollger/repos/fibertools-rs/benchmarks/2025-07-08-new-torch-version/bin/ft_0_17_1_pytorch_all_features "$@"