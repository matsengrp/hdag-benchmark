#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

ml Beast/2.6.3-GCC-10.2.0

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
ctreevcf=${ctreefasta_with_refseq}.vcf

beastdir=$resultsdir/beast
mkdir -p $beastdir

# TODO: Figure out how to initalize from the same start tree or at least figure out how to fix the ancestral sequence
#       Can we just outgroup using the ancestral sequence?
# starttree=$beastdir/starttree.nwk
# echo $dagdir
# echo "()ancestral;" > "$starttree"
# refseqfile=$dagdir/refseq.txt
# head -2 $ctreefasta_with_refseq | tail -n +2 > $refseqfile
# seedtree=$dagdir/seedtree.pb
# echo "===> create tree with UShER..."
# # $ctreevcf contains the ancestral sequence and all the other simulated sequences
# usher-sampled -v $ctreevcf -t $starttree -o $seedtree  --optimization_minutes=0

cd $beastdir



tree_file=$beastdir/beast-output.trees

# TODO: Need to uncomment this when we do a new beast run
# echo "===> Running beast..."
# # Create XML file
# conda activate beast-xml
# xml_file=default.xml
# beast2-xml.py \
# --fastaFile $ctreefasta_with_refseq > $xml_file \
# --chainLength 10000000000 \
# --treeLogEvery 1000
# # Run beast
# beast -overwrite $xml_file

# Going back to (what should be) the data directory
cd $currdir

# Extract support from trees
echo "===> Extracting supports..."
conda activate hdag-benchmark
python ../support_pipeline_scripts/cli.py save_supports -m "beast" -t $ctree -i $tree_file -o $beastdir/results.pkl
echo ""
echo ""

# python support_pipeline_scripts/cli.py save_supports -m "beast" -t "/home/whowards/hdag-benchmark/data/A.2.2/1/simulation/collapsed_simulated_tree.nwk" -i "/home/whowards/hdag-benchmark/data/A.2.2/1/results/beast/beast-output.trees" -o "/home/whowards/hdag-benchmark/data/A.2.2/1/results/beast/results_temp.pkl"