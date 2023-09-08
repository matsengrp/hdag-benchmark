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
sim_model=$4

echo $currdir
datadir=$currdir/data/sim_models

# cd into the data directory
cd $datadir

trialdir=$clade/$sim_model/$trial
simdir=$trialdir/"simulation"
resultsdir=$trialdir/"results"
mrbayesdir=$resultsdir/mrbayes
mkdir -p $mrbayesdir

ctree=$simdir/collapsed_simulated_tree.nwk
ctreefasta=${ctree}.fasta
ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta
ctreenexus=$mrbayesdir/ctree_with_refseq.nex
mrbayesfile=$mrbayesdir/run.mb

## ---- Comment this section out if you don't want to rerun mrbayes ---- ##

# # Convert fasta to nexus file for mrbayes
# seqmagick convert $ctreefasta_with_refseq $ctreenexus --alphabet dna
# mrbayesoutput=$mrbayesdir/mrbayes-output

# # Produce .mb file describing the mrbayes run (including input and output files)
# python $currdir/support_pipeline/scripts/python_replace.py \
#     $currdir/support_pipeline/scripts/mb_simulation_experiment/mb_files/run_unrest_hmut.mb \
#     $ctreenexus $mrbayesoutput > $mrbayesfile
# /fh/fast/matsen_e/whowards/MrBayes/src/mb -i $mrbayesfile

## --------------------------------------------------------------------- ##

# Although the mrbayes-output.trprobs file contains the deduplicated
# topologies, annotated with their posterior probabilities.
tree_file=$mrbayesdir/mrbayes-output.t


# Extract support from trees
rtree=$simdir/resolved_output.nwk
echo "===> Extracting supports..."
conda activate hdag-benchmark
python $currdir/support_pipeline/inference.py save_supports -m "mrbayes" \
-t $rtree \
-f $ctreefasta \
-i $tree_file \
-o $mrbayesdir/results.pkl
# --use_results

# python $currdir/support_pipeline/scripts/mb_simulation_experiment/mb_parsimony_distribution.py $sim_model $clade $trial


# python support_pipeline/inference.py save_supports -m "mrbayes" \
# -t /fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation/AY.108/0.01_200/1/simulation/resolved_output.nwk \
# -f /fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation/AY.108/0.01_200/1/simulation/collapsed_simulated_tree.nwk.fasta \
# -i /fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation/AY.108/0.01_200/1/results/mrbayes/mrbayes-output.t \
# -o /fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation/AY.108/0.01_200/1/results/mrbayes/results.pkl \
# --use_results

