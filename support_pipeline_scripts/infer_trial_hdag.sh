#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

# NOTE: This script assumes we are currently in the data directory
currdir=$1
clade=$2
trial=$3

trialdir=$currdir/$clade/$trial
simdir=$trialdir/"simulation"
resultsdir=$trialdir/"results"
mkdir -p $resultsdir

ctree=$simdir/collapsed_simulated_tree.nwk
ctreefasta=${ctree}.fasta
ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta
ctreevcf=${ctreefasta_with_refseq}.vcf

dagdir=$resultsdir/historydag
mkdir -p $dagdir
starttree=$dagdir/starttree.nwk
echo $dagdir
echo "()ancestral;" > "$starttree"
refseqfile=$dagdir/refseq.txt
head -2 $ctreefasta_with_refseq | tail -n +2 > $refseqfile
seedtree=$dagdir/seedtree.pb

mkdir -p $dagdir/opt_info
# cd $dagdir/opt_info

echo "===> create tree with UShER..."
# $ctreevcf contains the ancestral sequence and all the other simulated sequences
usher-sampled -v $ctreevcf -t $starttree -o $seedtree  --optimization_minutes=0 -d $dagdir/opt_info


echo "===> lusher optimizing..."
log_prefix=$dagdir/opt_info/optimization_log
optdag_final=$dagdir/final_opt_dag.pb


# NOTE this script changes directory to path/to/results/historydag
python ../support_pipeline_scripts/cli.py larch_usher -i $seedtree -r $refseqfile -c 2000 -o $dagdir -l $log_prefix
python ../support_pipeline_scripts/cli.py save_supports -m "hdag" -t $ctree -i $optdag_final -o $dagdir/results.pkl
echo ""
echo ""