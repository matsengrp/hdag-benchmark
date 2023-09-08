#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

# ml MrBayes/3.2.7a-foss-2021b

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
mrbayesdir=$resultsdir/mrbayes
mkdir -p $mrbayesdir

ctree=$simdir/collapsed_simulated_tree.nwk
ctreefasta=${ctree}.fasta
ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta
ctreenexus=$mrbayesdir/ctree_with_refseq.nex
mrbayesfile=$mrbayesdir/run.mb

# # Convert fasta to nexus file for mrbayes
# seqmagick convert $ctreefasta_with_refseq $ctreenexus --alphabet dna
# mrbayesoutput=$mrbayesdir/mrbayes-output

# # Produce .mb file describing the mrbayes run (including input and output files)
# python $currdir/support_pipeline/scripts/python_replace.py $currdir/run.mb $ctreenexus $mrbayesoutput > $mrbayesfile
# /fh/fast/matsen_e/whowards/MrBayes/src/mb -i $mrbayesfile

# Although the mrbayes-output.trprobs file contains the deduplicated
# topologies, annotated with their posterior probabilities.
tree_file=$mrbayesdir/mrbayes-output.t

# Going back to (what should be) the data directory
cd $datadir

# Extract support from trees
echo "===> Extracting supports..."
python $currdir/support_pipeline/inference.py save_supports -m "mrbayes" -t $ctree -i $tree_file -o $mrbayesdir/results.pkl
echo ""
echo ""

