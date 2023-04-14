#!/bin/bash

# Aggregates result.pkl files for all trials. Assumes directory structure created by simulation
# and inference scripts
# Usage is aggregate_results <method> (where <method> is either beast or historydag)

# IMPORTANT: This script assumes it is being run from hdag-benchmark directory

set -eu
method=$1
num_res=10
num_sim=10
let "num_trials = $num_sim * $num_res"

echo ""
echo "=> Aggregating results..."


eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python


for clade in $(cat clades.txt); do
    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    cladedir=data/$clade

    for trial in $(seq $num_trials); do
        echo "$clade: $trial"
        outdir=$cladedir/$trial/figures/$method
        mkdir -p $outdir

        results=$cladedir/$trial/results/$method/results.pkl
        # python support_pipeline_scripts/cli.py agg -i $results -o $outdir -c $clade -w 0.2 -m $method
        
        # results=$cladedir/$trial/results/$method/strat_dict_pars_weight.pkl
        # python support_pipeline_scripts/cli.py agg_pars_weights -i $results -o $outdir -c $clade -w 0.2 -m $method
        # python support_pipeline_scripts/cli.py bin_pars_weights -i $results -o $outdir -c $clade -b 0.05 -m $method
        dag_path=$cladedir/$trial/results/$method/final_opt_dag.pb
        python support_pipeline_scripts/cli.py cumm_pars_weight -i $dag_path -o $outdir -p 0.1

    done

    outdir=$cladedir/figures/$method
    mkdir -p $outdir
    # python support_pipeline_scripts/cli.py clade_results -n $num_trials -c $cladedir -o $outdir -m $method

    python support_pipeline_scripts/cli.py pars_weight_clade_results -n $num_trials -c $cladedir -o $outdir -m $method
done

# python support_pipeline_scripts/cli.py agg \
# -i "/home/whowards/hdag-benchmark/data/A.2.2/1/results/beast/results.pkl" \
# -o "/home/whowards/hdag-benchmark/data/A.2.2/1/figures/beast" \
# -c A.2.2 \
# -w 0.2 \
# -m beast
#