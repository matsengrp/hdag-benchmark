#!/bin/bash

# Uses ...

num_res=2
num_sim=2
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
        results=$cladedir/$trial/results/historydag/results.pkl
        outdir=$cladedir/$trial/figures
        mkdir -p $outdir

        python ~/hdag-benchmark/support_pipeline_scripts/cli.py agg -i $results -o $outdir
    done

    # TODO: Aggregate all results into a single plot
    #       - This could be done in a couple of ways
    #       (1) Append all results list together and do the sliding window on a much larger list
    #       (2) Do sliding window on four separate lines, plot all 4 (or take the average somehow)
    mkdir -p $cladedir/figures
    python ~/hdag-benchmark/support_pipeline_scripts/cli.py clade_results -n $num_trials -c $cladedir
done
