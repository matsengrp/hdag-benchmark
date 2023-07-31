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

currdir=$PWD
datadir=data/hypermutation
script_dir=$currdir"/support_pipeline/scripts/hypermutation_experiment"
num_trials=3

for clade in $(cat $script_dir/clades.txt); do
    # if value of $clade starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade

    for hmut in $(cat $script_dir/hypermutation.txt); do

        [[ $hmut =~ ^#.* ]] && continue

        echo "   $hmut"

        cladedir=$datadir/$clade/$hmut
        outdir=$cladedir/figures/$method
        mkdir -p $outdir

        python support_pipeline/plotting.py clade_results \
        -n $num_trials \
        -c $cladedir \
        -m $method \
        -r results.pkl \
        -o $outdir/CA_support.png
    
    done

done

# python support_pipeline/plotting.py clade_results \
# -n 3 \
# -c "/fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation/AY.108/0.01_200" \
# -m mrbayes \
# -r "results.pkl" \
# -o "/fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation/AY.108/0.01_200/figures/mrbayes/CA_support.png"
