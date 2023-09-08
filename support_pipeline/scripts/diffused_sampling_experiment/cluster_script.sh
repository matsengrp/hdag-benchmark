#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"

conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
export QT_QPA_PLATFORM=offscreen; export XDG_RUNTIME_DIR=/tmp/runtime-runner; export MPLBACKEND=agg; ml libGLU/9.0.1-GCCcore-10.2.0

# NOTE: This file must have executable permissions to be used with cluster driver

script_name=$1
clade=$2
trial=$3

echo $1
echo $2
echo $3

python support_pipeline/scripts/diffused_sampling_experiment/$script_name $clade $trial