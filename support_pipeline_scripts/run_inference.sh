#!/bin/bash

# Uses usher (and other methods later) to create a starting tree on the simulated data, then uses larch-usher to optimize
# This script is intended to be run after `run_simulations.sh`
# Should output file containing node supports for every simulated dataset

# Paramaters that determine how many trees to simulate for each clade
num_res=2
num_sim=2
let "num_trials = $num_sim * $num_res"
# TODO: makes these parameters to pipeline script

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
    # TODO: send a cluster job off that does all this stuff?
    echo $clade
    cladedir=$datadir/$clade

    for trial in $(seq $num_trials); do
        trialdir=$cladedir/$trial
        simdir=$trialdir/"simulation"
        resultsdir=$trialdir/"results"
        mkdir -p $resultsdir

        ctree=$simdir/collapsed_simulated_tree.nwk
        ctreefasta=${ctree}.fasta
        ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta
        ctreevcf=${ctreefasta_with_refseq}.vcf

        dagdir=$resultsdir/historydag
        mkdir -p $dagdir
        starttree=$dagdir/starttree.nwk
        echo $dagdir
        echo "()ancestral;" > "$starttree"
        refseqfile=$dagdir/refseq.txt
        head -2 $ctreefasta_with_refseq | tail -n +2 > $refseqfile
        seedtree=$dagdir/seedtree.pb

        mkdir -p $dagdir/opt_info
        cd $dagdir/opt_info

        echo "===> create tree with UShER..."
        # $ctreevcf contains the ancestral sequence and all the other simulated sequences
        usher-sampled -v $ctreevcf -t $starttree -o $seedtree  --optimization_minutes=0

        echo "===> lusher optimizing..."
        optimized_dag=$dagdir/opt_dag.pb
        $larch_usher_exec -i $seedtree -r $refseqfile -c 100 -o $optimized_dag

        python ~/hdag-benchmark/support_pipeline_scripts/cli.py save_supports -m "hdag" -t $ctree -i $optimized_dag -o $dagdir/results.pkl
    done
done
