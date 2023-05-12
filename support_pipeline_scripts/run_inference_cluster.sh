#!/bin/bash

# Uses usher (and other methods later) to create a starting tree on the simulated data, then uses larch-usher to optimize
# This script is intended to be run after `run_simulations.sh`
# Should output file containing node supports for every simulated dataset

# Paramaters that determine how many trees to simulate for each clade
num_res=10
num_sim=10
method="historydag"
let "num_trials = $num_sim * $num_res"

echo ""
echo "=> Starting inference..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python


currdir=$PWD
datadir=data
cd $datadir

for clade in $(cat ../clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade
    cladedir=$datadir/$clade

    for trial in $(seq $num_trials); do
        logfile=$clade/$trial/results/inference.log
        echo $logfile $method

        mkdir -p $clade/$trial/results
        
        sbatch -c 1 -J "$trial|$clade|inference" -o $logfile -e $logfile.err \
        $currdir/support_pipeline_scripts/infer_trial_$method.sh $currdir $clade $trial

    done
done
