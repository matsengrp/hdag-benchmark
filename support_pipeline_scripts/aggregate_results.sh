#!/bin/bash

# Aggregates result.pkl files for all trials. Assumes directory structure created by simulation
# and inference scripts
# Usage is aggregate_results <method> (where <method> is either beast or historydag)

method=$1
num_res=10
num_sim=10
let "num_trials = $num_sim * $num_res"

echo ""
echo "=> Aggregating results..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

# 
cd ~/hdag-benchmark/data

for clade in $(cat ../clades.txt); do
    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    cladedir=$clade

    for trial in $(seq $num_trials); do
        echo "$clade: $trial"
        results=$cladedir/$trial/results/$method/results.pkl
        outdir=$cladedir/$trial/figures/$method
        mkdir -p $outdir

        python ~/hdag-benchmark/support_pipeline_scripts/cli.py agg -i $results -o $outdir -c $clade -w 0.2 -m $method
    done

    outdir=$cladedir/figures/$method
    mkdir -p $outdir
    python ~/hdag-benchmark/support_pipeline_scripts/cli.py clade_results -n $num_trials -c $cladedir -o $outdir -m $method
done
