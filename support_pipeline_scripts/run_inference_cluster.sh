#!/bin/bash

# Uses usher (and other methods later) to create a starting tree on the simulated data, then uses larch-usher to optimize
# This script is intended to be run after `run_simulations.sh`
# Should output file containing node supports for every simulated dataset

# Paramaters that determine how many trees to simulate for each clade
num_res=10
num_sim=10
method="beast"
let "num_trials = $num_sim * $num_res"

echo ""
echo "=> Starting inference..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

# Path to larch-usher executable (needs to be pre-built)
larch_usher_exec=/home/whowards/larch/larch/build/larch-usher


datadir=data
cd $datadir

for clade in $(cat ../clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade
    cladedir=$datadir/$clade

    for trial in $(seq $num_trials); do
        logfile=$clade/$trial/results/inference.log
        echo $logfile

        mkdir -p $clade/$trial/results
        
        sbatch -c 4 -J "$trial|$clade|inference" -o $logfile -e $logfile.err \
        ./../support_pipeline_scripts/infer_trial_$method.sh ~/hdag-benchmark/data $clade $trial
    done
done