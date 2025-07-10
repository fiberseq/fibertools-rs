#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate pytorch2.2
export LIBTORCH_USE_PYTORCH=1
export DYLD_LIBRARY_PATH=/Users/mrvollger/miniconda3/envs/pytorch2.2/lib/python3.9/site-packages/torch/lib/:$DYLD_LIBRARY_PATH
./ft_burn_0_12_pytorch "$@"

