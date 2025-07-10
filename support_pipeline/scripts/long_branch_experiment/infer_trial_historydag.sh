#!/bin/bash

source /n/fs/ragr-research/users/wh8114/anaconda3/bin/activate hdag-benchmark

conda activate hdag-benchmark

# NOTE: This script assumes we are currently in the data/sim_models directory
cd /n/fs/ragr-research/users/wh8114/prev/hdag-benchmark/data/sim_models
currdir=$1
trialdir=$2

simdir=$trialdir/"simulation"
resultsdir=$trialdir/"results"
echo "Making $resultsdir"
mkdir -p $resultsdir

ctree=$simdir/collapsed_simulated_tree.nwk
ctreefasta=${ctree}.fasta
ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta
ctreevcf=${ctreefasta_with_refseq}.vcf

dagdir=$resultsdir/historydag
echo "Making $dagdir/opt_info"
mkdir -p $dagdir/opt_info

starttree=$dagdir/starttree.nwk
echo "()ancestral;" > "$starttree"
refseqfile=$dagdir/refseq.txt
head -2 $ctreefasta_with_refseq | tail -n +2 > $refseqfile
seedtree=$dagdir/seedtree.pb_tree

optdag_final=final_opt_dag_trimmed

echo "===> create tree with UShER..."
# $ctreevcf contains the ancestral sequence and all the other simulated sequences
usher-sampled -v $ctreevcf -t $starttree -o $seedtree --optimization_minutes=0 -d $dagdir/opt_info

# conda activate larch-env
source /n/fs/ragr-research/users/wh8114/anaconda3/bin/activate larch

echo "===> lusher optimizing..."
log_prefix=$dagdir/opt_info/optimization_log
mkdir -p $log_prefix

/n/fs/ragr-research/users/wh8114/prev/larch/build/bin/larch-usher \
-i $seedtree \
-o $dagdir/final_opt_dag_trimmed.pb \
-r $refseqfile \
--move-coeff-nodes 1 \
--move-coeff-pscore 5 \
--sample-any-tree \
--trim \
-l $log_prefix \
--inter-save 10 \
-c 1000 \
--input-format tree-pb \
--output-format dag-pb