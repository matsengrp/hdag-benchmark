#!/bin/bash

# Uses usher (and other methods later) to create a starting tree on the simulated data, then uses larch-usher to optimize
# This script is intended to be run after `run_simulations.sh`
# Should output file containing node supports for every simulated dataset

# Paramaters that determine how many trees to simulate for each clade
num_res=1
num_sim=1
let "num_trials = $num_sim * $num_res"

echo ""
echo "=> Starting inference..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

# Path to larch-usher executable (needs to be pre-built)
larch_usher_exec=/home/whowards/larch/larch/build/larch-usher

# mkdir -p ~/hdag-benchmark/opt_logs
# cd ~/hdag-benchmark/opt_logs

datadir=~/hdag-benchmark/data
for clade in $(cat ../clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade
    cladedir=$datadir/$clade

    for trial in $(seq $num_trials); do
        echo $cladedir/$trial/results/inference.log
        sbatch -c 1 -J "$clade|inf|$trial" -o $cladedir/$trial/results/inference.log \
        ./../support_pipeline_scripts/infer_trial.sh ~/hdag-benchmark/data $clade $trial
    done
done

# python ~/hdag-benchmark/support_pipeline_scripts/cli.py save_supports -m "hdag" \
# -t /home/whowards/hdag-benchmark/data/A.2.5/1/simulation/collapsed_simulated_tree.nwk \
# -i /home/whowards/hdag-benchmark/data/A.2.5/1/results/historydag/opt_dag.pb \
# -o /home/whowards/hdag-benchmark/data/A.2.5/1/results/historydag/results.pkl