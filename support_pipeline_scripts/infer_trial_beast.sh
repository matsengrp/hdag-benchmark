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

echo "===> Running beast..."

# TODO: Should we be using the fasta with ancestral sequence???
# TODO: Set beast xml params correctly
#       - You can use template file that has correct model settings
#       - What chain length do you want?

tree_file=$beastdir/beast-output.trees

# Create XML file
conda activate beast-xml
xml_file=default.xml
beast2-xml.py \
--fastaFile $ctreefasta > $xml_file \
--chainLength 100000000 \
--treeLogEvery 1000

# Run beast
beast -overwrite $xml_file

# Extract support from trees
# TODO: Change this to not have the 2 at the end
conda activate hdag-benchmark
python ~/hdag-benchmark/support_pipeline_scripts/cli.py save_supports -m "beast" -t $ctree -i $tree_file -o $beastdir/results.pkl
echo ""
echo ""