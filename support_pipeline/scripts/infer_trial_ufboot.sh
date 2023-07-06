#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

currdir=$1
clade=$2
trial=$3

echo $currdir
datadir=$currdir/data

# cd into the data directory
cd $datadir

trialdir=$datadir/$clade/$trial
simdir=$trialdir/"simulation"
resultsdir=$trialdir/"results"
ufbootdir=$resultsdir/ufboot
mkdir -p $ufbootdir

ctree=$simdir/collapsed_simulated_tree.nwk
ctreefasta=${ctree}.fasta
ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta

# cd $ufbootdir
# iqtree -s $ctreefasta_with_refseq -B 1000 --bnni --prefix output

out_file=$ufbootdir"/output.splits.nex"

# Going back to (what should be) the data directory
cd $datadir

# Extract support from trees
echo "===> Extracting supports..."
conda activate hdag-benchmark
python $currdir/support_pipeline/inference.py save_supports -m "ufboot" -t $ctree -i $out_file -o $ufbootdir/results.pkl
echo ""
echo ""

