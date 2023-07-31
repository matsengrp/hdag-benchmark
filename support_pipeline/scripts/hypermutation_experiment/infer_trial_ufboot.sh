#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

currdir=$1
clade=$2
trial=$3
hmut=$4

echo $currdir
datadir=$currdir/data/hypermutation

# cd into the data directory
cd $datadir

trialdir=$datadir/$clade/$hmut/$trial
simdir=$trialdir/"simulation"
resultsdir=$trialdir/"results"
ufbootdir=$resultsdir/ufboot
mkdir -p $ufbootdir

ctree=$simdir/collapsed_simulated_tree.nwk
ctreefasta=${ctree}.fasta
ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta

cd $ufbootdir

# Uncomment this line to rerun inference
# iqtree -s $ctreefasta_with_refseq -B 1000 --bnni --prefix output

out_file=$ufbootdir"/output.splits.nex"

# Going back to (what should be) the data directory
cd $datadir

# Extract support from trees
echo "===> Extracting supports..."
conda activate hdag-benchmark
rtree=$simdir/resolved_output.nwk
python $currdir/support_pipeline/inference.py save_supports -m "ufboot" -t $rtree -i $out_file -f $ctreefasta -o $ufbootdir/results.pkl --use_results
echo ""
echo ""
