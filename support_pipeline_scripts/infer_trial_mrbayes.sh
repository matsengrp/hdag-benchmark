#!/bin/bash
# NOTE: This file must have executable permissions to be used with simulation pipeline

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python


# ml MrBayes/3.2.7a-gompi-2021b

currdir=$1
clade=$2
trial=$3
echo $trial

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

echo converting fasta to nexus
# Convert fasta to nexus file for mrbayes
seqmagick convert $ctreefasta_with_refseq $ctreenexus --alphabet dna
mrbayesoutput=$mrbayesdir/mrbayes-output

echo getting true tree
# get the resolved true tree newick with branch lengths:
resolvedtruetree=$mrbayesdir/resolved_true_tree.nwk
python $currdir/hdb/cli.py resolve-multifurcations -i $ctree -o $resolvedtruetree --add-ancestral

scaledresolvedtree=$mrbayesdir/scaled_resolved_tree.nwk
python $currdir/support_pipeline_scripts/cli.py scale_branch_lengths -i $resolvedtruetree -s 0.0000345 > $scaledresolvedtree

echo building mrbayes file
# Produce .mb file describing the mrbayes run (including input and output files)
python $currdir/support_pipeline_scripts/python_replace.py $currdir/run.mb $ctreenexus $mrbayesoutput > $mrbayesfile

/fh/fast/matsen_e/whowards/MrBayes/src/mb -i $mrbayesfile

# Although the mrbayes-output.trprobs file contains the deduplicated
# topologies, annotated with their posterior probabilities.
tree_file=$mrbayesdir/mrbayes-output.t

# Going back to (what should be) the data directory
cd $datadir

# Extract support from trees
echo "===> Extracting supports..."
conda activate hdag-benchmark
python $currdir/support_pipeline_scripts/cli.py save_supports -m "mrbayes" -t $ctree -i $tree_file -o $mrbayesdir/results.pkl
echo ""
echo ""

# python /fh/fast/matsen_e/whowards/hdag-benchmark/support_pipeline_scripts/cli.py save_supports \
# -m "mrbayes" \
# -t /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.108/16/simulation/collapsed_simulated_tree.nwk \
# -i /fh/fast/matsen_e/wdumm/mrbayes-test/test_output.t \
# -o /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.108/16/results/mrbayes/results.pkl

