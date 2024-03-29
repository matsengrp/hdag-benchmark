#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"

conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
export QT_QPA_PLATFORM=offscreen; export XDG_RUNTIME_DIR=/tmp/runtime-runner; export MPLBACKEND=agg; ml libGLU/9.0.1-GCCcore-10.2.0

# NOTE: This file must have executable permissions to be used with cluster driver

clade=$1
trial=$2
rf_method=$3

echo $1
echo $2
echo $3

python support_pipeline/scripts/rf_dist_scaling_experiment/rf_scaling.py $clade $trial $rf_method

# E.g.
# bash support_pipeline/scripts/rf_dist_scaling_experiment/rf_scaling.sh \
# AY.34.2 \
# 1 \
# min_trim
