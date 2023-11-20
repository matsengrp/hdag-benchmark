#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

# NOTE: This script assumes we are currently in the data directory
currdir=$1
clade=$2
trial=$3
sim_model=$4

python /fh/fast/matsen_e/whowards/hdag-benchmark/support_pipeline/scripts/mb_simulation_experiment/mb_parsimony_distribution.py\
 $sim_model $clade $trial