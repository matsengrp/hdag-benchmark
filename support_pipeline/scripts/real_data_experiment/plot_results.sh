#!/bin/bash

# NOTE: These intructions are for the original version of this script and have not been updated
# Aggregates result.pkl files for all trials. Assumes directory structure created by simulation
# and inference scripts
# Usage is bash plot_results <method> (where <method> is either mrbayes or historydag)

# IMPORTANT: This script assumes it is being run from hdag-benchmark directory

set -eu
method=$1


echo ""
echo "=> Plotting results..."


eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

for clade in $(cat clades.txt); do
    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    cladedir=data/$clade"_"

    echo "$clade"
    outdir=$cladedir/figures/$method
    mkdir -p $outdir

    results=$cladedir/results/$method/results_adj_first.pkl
    echo $results
    python support_pipeline/plotting.py coverage_trial_plot -i $results -o $outdir -c $clade -w 0.2 -m $method

done
