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

datadir=$currdir/data
trialdir=$datadir/$clade/$trial
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
optdag_final=final_opt_dag

# TODO: Uncomment this when you want to search for new trees ---v
# $ctreevcf contains the ancestral sequence and all the other simulated sequences
echo "===> create tree with UShER..."
usher-sampled -v $ctreevcf -t $starttree -o $seedtree  --optimization_minutes=0 -d $dagdir/opt_info

echo "===> lusher optimizing..."
log_prefix=$dagdir/opt_info/optimization_log
python ../support_pipeline/inference.py larch_usher -i $seedtree -r $refseqfile -c 4000 -o $dagdir -l $log_prefix -f $optdag_final

# TODO: Rename 'hdag-...' -> 'historydag-...' and change save_supports method accordingly
# NOTE: Saves single results list for MP trimmed node support
# python ../support_pipeline/inference.py save_supports -m "hdag-adj" -t $ctree -i $dagdir/$optdag_final.pb -o $dagdir/adjusted_support_results.pkl
# TODO: Change back to optdag_final
python ../support_pipeline/inference.py save_supports -m "hdag-inf" -t $ctree -i $dagdir/$optdag_final.pb -o $dagdir/results.pkl

# NOTE: Tests different trimming strategies for uniform node support.
# python ../support_pipeline/inference.py trim_thresholds -t $ctree -i $optdag_final -o $dagdir/strat_dict.pkl

# NOTE: Tests different weighting schemes for generalized node support.
# python ../support_pipeline/inference.py test_pars_weights -t $ctree -i $dagdir/$optdag_final.pb -o $dagdir/strat_dict_node_weight.pkl
echo ""
echo ""
