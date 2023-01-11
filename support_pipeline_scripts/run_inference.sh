#!/bin/bash

# Uses usher (and other methods later) to create a starting tree on the simulated data, then uses larch-usher to optimize
# This script is intended to be run after `run_simulations.sh`
# Should output file containing node supports for every simulated dataset

echo ""
echo "=> Starting inference..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

# Path to larch-usher executable (needs to be pre-built)
larch_usher_exec=/home/whowards/larch/larch/build/larch-usher

cd ~/hdag-benchmark/data

ignoredir=_ignore
mkdir -p $ignoredir
for clade in $(cat clades.txt); do
    # TODO: send a cluster job off that does all this stuff?
    echo $clade
    cladedir=$ignoredir/$clade

    ctree=$cladedir/collapsed_simulated_tree.nwk
    ctreefasta=${ctree}.fasta
    ctreefasta_with_refseq=$cladedir/ctree_with_refseq.fasta
    ctreevcf=${ctreefasta_with_refseq}.vcf

    dagdir=$cladedir/historydag
    mkdir -p $dagdir
    starttree=$dagdir/starttree.nwk
    echo "()ancestral;" > $starttree
    refseqfile=$dagdir/refseq.txt
    head -2 $ctreefasta_with_refseq | tail -n +2 > $refseqfile
    seedtree=$dagdir/seedtree.pb

    echo "===> create tree with UShER..."
    # $ctreevcf contains the ancestral sequence and all the other simulated sequences
    usher-sampled -v $ctreevcf -t $starttree -o $seedtree  --optimization_minutes=0

    echo "===> lusher optimizing..."
    optimized_dag=$dagdir/opt_dag.pb
    $larch_usher_exec -i $seedtree -r $refseqfile -c 1000 -o $optimized_dag

    # TODO: Open DAG in Python and run node-support then save as file.
    python ~/hdag-benchmark/support_pipeline_scripts/cli.py save_supports -m "hdag" -t $ctree -i $optimized_dag -o $cladedir/results.pkl
done
