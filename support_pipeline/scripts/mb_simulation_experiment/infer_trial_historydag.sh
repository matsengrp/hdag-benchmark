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
sim_model=$4

datadir=$currdir/data/sim_models
trialdir=$datadir/$clade/$sim_model/$trial
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
optdag_final=final_opt_dag_trimmed

# TODO: Uncomment this when you want to search for new trees ---v
# $ctreevcf contains the ancestral sequence and all the other simulated sequences
# echo "===> create tree with UShER..."
# usher-sampled -v $ctreevcf -t $starttree -o $seedtree --optimization_minutes=0 -d $dagdir/opt_info

# echo "===> lusher optimizing..."
# log_prefix=$dagdir/opt_info/optimization_log
# python $currdir/support_pipeline/inference.py larch_usher -i $seedtree -r $refseqfile -c 5000 -o $dagdir -l $log_prefix -f $optdag_final

# python $currdir/support_pipeline/inference.py save_supports -m "hdag" -t $ctree -i $dagdir/$optdag_final.pb -o $dagdir/results.pkl -f $ctreefasta # --use_results

# Sample from diffused dag distribution and save in input file
diff_dagdir=$resultsdir/diff_historydag_mut_rates
mkdir -p $diff_dagdir
python $currdir/support_pipeline/inference.py diffused_hdag_samples -o $diff_dagdir/diff_sample.trees \
-i $dagdir/$optdag_final.pb -n 100000 -f $ctreefasta -m $simdir/sars-cov-2_simulation_output.info
python $currdir/support_pipeline/inference.py save_supports -m "hdag-diff" -t $ctree -i $diff_dagdir/diff_sample.trees -o $diff_dagdir/results.pkl -f $ctreefasta # --use_results
