#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

currdir=$1
clade=$2
trial=1

echo $currdir
datadir=$currdir/data

# cd into the data directory
cd $datadir

cladedir=$datadir/$clade"_"
trialdir=$cladedir/trials/$trial
resultsdir=$cladedir/"results"
ufbootdir=$resultsdir/ufboot
mkdir -p $ufbootdir

fasta_file=$cladedir/reconstructed_seqs.fasta

## ---- Comment this section out if you don't want to rerun ufboot ---- ##

# cd $ufbootdir
# iqtree -s $fasta_file -B 1000 --bnni --prefix output
# cd $datadir

## --------------------------------------------------------------------- ##

out_file=$ufbootdir"/output.splits.nex"

python $currdir/support_pipeline/inference.py save_supports \
-m "ufboot" \
-t $trialdir/sampled_tree.nwk \
-i $out_file \
-o $ufbootdir/results.pkl