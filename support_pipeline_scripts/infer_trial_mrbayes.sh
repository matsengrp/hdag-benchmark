#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

ml MrBayes/3.2.7a-foss-2021b

currdir=$1
clade=$2
trial=$3

echo $currdir

# cd into the data directory
cd $currdir

trialdir=$currdir/$clade/$trial
simdir=$trialdir/"simulation"
resultsdir=$trialdir/"results"
mkdir -p $resultsdir

ctree=$simdir/collapsed_simulated_tree.nwk
ctreefasta=${ctree}.fasta
ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta
ctreenexus=$simdir/ctree_with_refseq.nex

# Convert fasta to nexus file for mrbayes
seqmagick convert $ctreefasta_with_refseq $ctreenexus --alphabet dna
mrbayesdir=$resultsdir/mrbayes
mkdir -p $mrbayesdir
mrbayesfile=$mrbayesdir/run.mb
mrbayesoutput="mrbayes-output"
# Produce .mb file describing the mrbayes run (including input and output files)
python support_pipeline_scripts/python_replace.py ../../run.mb $ctreefasta_with_refseq $mrbayesoutput > $mrbayesfile

mb -i $mrbayesfile

# Although the mrbayes-output.trprobs file contains the deduplicated
# topologies, annotated with their posterior probabilities.
tree_file=$mrbayesdir/mrbayes-output.trees

# Going back to (what should be) the data directory
cd $currdir

# Extract support from trees
echo "===> Extracting supports..."
conda activate hdag-benchmark
python ../support_pipeline_scripts/cli.py save_supports -m "mrbayes" -t $ctree -i $tree_file -o $mrbayesdir/results.pkl
echo ""
echo ""

# python support_pipeline_scripts/cli.py save_supports -m "beast" -t "/home/whowards/hdag-benchmark/data/A.2.2/1/simulation/collapsed_simulated_tree.nwk" -i "/home/whowards/hdag-benchmark/data/A.2.2/1/results/beast/beast-output.trees" -o "/home/whowards/hdag-benchmark/data/A.2.2/1/results/beast/results_temp.pkl"