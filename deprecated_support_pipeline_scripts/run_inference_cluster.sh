#!/bin/bash

# Uses usher (and other methods later) to create a starting tree on the simulated data, then uses larch-usher to optimize
# This script is intended to be run after `run_simulations.sh`
# Should output file containing node supports for every simulated dataset

# Paramaters that determine how many trees to simulate for each clade
num_res=10
num_sim=10
let "num_trials = $num_sim * $num_res"

num_cores=4
trial_file="pars_div_trials"
method=$1


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

    for trial in $(cat ../$trial_file.txt); do
        logfile=$clade/$trial/results/inference_$method.log
        echo $logfile $method

        mkdir -p $clade/$trial/results
        
        sbatch -t 6-0 -c $num_cores -J "$trial|$clade|inference" -o $logfile -e $logfile.err \
        $currdir/support_pipeline_scripts/infer_trial_$method.sh $currdir $clade $trial

    done
done