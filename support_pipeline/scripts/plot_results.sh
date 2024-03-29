#!/bin/bash

# Aggregates result.pkl files for all trials. Assumes directory structure created by simulation
# and inference scripts
# Usage is bash plot_results <method> (where <method> is either mrbayes or historydag)

# IMPORTANT: This script assumes it is being run from hdag-benchmark directory

set -eu
method=$1
sim_model="gamma_10_hmut_50"

num_res=25
num_sim=1
let "num_trials = $num_sim * $num_res"

echo ""
echo "=> Plotting results..."


eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

for clade in $(cat /fh/fast/matsen_e/whowards/hdag-benchmark/support_pipeline/scripts/diffused_sampling_experiment/clades.txt); do
    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    cladedir=data/sim_models/$clade/gamma_10_hmut_50

    # TODO: Make a script for plotting clade results separate from trial results
    
    for trial in $(seq 25); do
        echo "$clade: $trial"
        outdir=$cladedir/$trial/figures/$method
        mkdir -p $outdir

        # results=$cladedir/$trial/results/$method/results.pkl
        # python support_pipeline/plotting.py coverage_trial_plot -i $results -o $outdir -c $clade -w 0.2 -m $method
        
        # results=$cladedir/$trial/results/$method/strat_dict_pars_weight.pkl
        # python support_pipeline/plotting.py agg_pars_weights -i $results -o $outdir -c $clade -w 0.2 -m $method
        # python support_pipeline/plotting.py bin_pars_weights -i $results -o $outdir -c $clade -b 0.05 -m $method

        # dag_path=$cladedir/$trial/results/$method/final_opt_dag.pb
        # python support_pipeline/plotting.py cumul_pars_weight -i $dag_path -o $outdir -p -0.1

    done

    outdir=$cladedir/figures/$method
    mkdir -p $outdir

    # Convert space seperated list of trials to comma separated
    # echo $t | tr " " "," > script_temp.txt
    # t_comma=$(cat script_temp.txt)
    # rm script_temp.txt

    python support_pipeline/plotting.py clade_results \
    -n $num_trials \
    -c $cladedir \
    -m $method \
    -r results.pkl \
    -o $outdir/CA_support.png
    # -t $t_comma

    # NOTE: This code generates estimated supports for adjustment
    # python support_pipeline_scripts/cli.py clade_results \
    # -n $num_trials \
    # -c $cladedir \
    # -m $method \
    # -r adjusted_support_results.pkl \
    # -o $outdir/CA_adjusted_support.png
    # -r results.pkl \
    # -o $outdir/CA_support.png

    # NOTE: This code randomly scales some estimated supports
    # sample_size=0.35
    # python support_pipeline_scripts/cli.py clade_results_random_scaling \
    # -n $num_trials \
    # -c $cladedir \
    # -m $method \
    # -r results.pkl \
    # -s $sample_size \
    # -o $outdir/CA_support_random_scaling_s=$sample_size-r=0.5-1.png

    # python support_pipeline_scripts/cli.py pars_weight_clade_results -n $num_trials -c $cladedir -o $outdir -m $method
done
