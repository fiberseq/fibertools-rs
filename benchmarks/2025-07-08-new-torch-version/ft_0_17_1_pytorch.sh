#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate pytorch2.6
export LIBTORCH_USE_PYTORCH=1
export DYLD_LIBRARY_PATH=/Users/mrvollger/miniconda3/envs/pytorch2.6/lib/python3.13/site-packages/torch/lib/:$DYLD_LIBRARY_PATH

./ft_burn_0_17_1_pytorch "$@"

