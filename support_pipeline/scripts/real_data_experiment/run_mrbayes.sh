#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

# ml MrBayes/3.2.7a-foss-2021b

currdir=$1
clade=$2
trial=1

echo $currdir

# TODO: Change this to get out of old_data soon
datadir=$currdir/data/real_data

# cd into the data directory
cd $datadir

cladedir=$datadir/$clade
resultsdir=$cladedir/"results"
mrbayesdir=$resultsdir/mrbayes
mkdir -p $mrbayesdir

fasta_file=$cladedir/reconstructed_seqs.fasta
nexus_file=$mrbayesdir/reconstructed_seqs.nex
mrbayesfile=$mrbayesdir/run.mb

## ---- Comment this section out if you don't want to rerun mrbayes ---- ##

# # Convert fasta to nexus file for mrbayes
# seqmagick convert $fasta_file $nexus_file --alphabet dna
# mrbayesoutput=$mrbayesdir/mrbayes-output

# # Produce .mb file describing the mrbayes run (including input and output files)
# python $currdir/support_pipeline/scripts/python_replace.py $currdir/run.mb $nexus_file $mrbayesoutput > $mrbayesfile
# /fh/fast/matsen_e/whowards/MrBayes/src/mb -i $mrbayesfile

## --------------------------------------------------------------------- ##

# Prepare directories and move files around
trialdir=$cladedir/trials/$trial
mkdir -p $trialdir
cp $cladedir/reconstructed_seqs.fasta $trialdir/sampled_tree.nwk.fasta

# NOTE: Only need to do this once
# python $currdir/support_pipeline/scripts/real_data_experiment/sample_from_mb.py $mrbayesdir/mrbayes-output.t $trialdir/sampled_tree.nwk

# python $currdir/support_pipeline/inference.py save_supports \
# -m "mrbayes" \
# -t $trialdir/sampled_tree.nwk \
# -i $mrbayesdir/mrbayes-output.t \
# -o $mrbayesdir/results.pkl

python $currdir/support_pipeline/scripts/sandbox/mb_parsimony_distribution.py $clade
