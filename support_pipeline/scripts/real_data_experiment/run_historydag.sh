#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE:
# - This file must have executable permissions to be used with simulation pipeline
# - This script assumes we are currently in the directory above `data`

currdir=$1
clade=$2
trial=1

datadir=$currdir/data
cladedir=$datadir/$clade"_"
trialdir=$cladedir/trials/$trial
resultsdir=$cladedir/"results"
mkdir -p $resultsdir

fasta_file=$cladedir/reconstructed_seqs.fasta
vcf_file=${fasta_file}.vcf
faToVcf $fasta_file $vcf_file

dagdir=$resultsdir/historydag
mkdir -p $dagdir
starttree=$dagdir/starttree.nwk
echo $dagdir
echo "()ancestral;" > "$starttree"
refseqfile=$dagdir/refseq.txt
head -2 $fasta_file | tail -n +2 > $refseqfile
seedtree=$dagdir/seedtree.pb

mkdir -p $dagdir/opt_info
optdag_final=final_opt_dag

## ---- Comment this section out if you don't want to rerun hdag ---- ##

# # $ctreevcf contains the ancestral sequence and all the other simulated sequences
# echo "===> create tree with UShER..."
# usher-sampled -v $vcf_file -t $starttree -o $seedtree  --optimization_minutes=0 -d $dagdir/opt_info

# echo "===> lusher optimizing..."
# log_prefix=$dagdir/opt_info/optimization_log
# python ../support_pipeline/inference.py larch_usher -i $seedtree -r $refseqfile -c 1000 -o $dagdir -l $log_prefix -f $optdag_final

## ------------------------------------------------------------------ ##

cp $dagdir/results.pkl $dagdir/results_copy.pkl

python $currdir/support_pipeline/inference.py save_supports \
-m "hdag" \
-t $trialdir/sampled_tree.nwk \
-i $dagdir/final_opt_dag.pb \
-o $dagdir/results.pkl \
--use_results
# NOTE: Remove the `--use_results` flag if you want to recompute the results list from scratch


# python support_pipeline/inference.py save_supports \
# -m "hdag-inf" \
# -t /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.108/25/simulation/collapsed_simulated_tree.nwk \
# -i /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.108/25/results/historydag/final_opt_dag.pb \
# -o /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.108/25/results/historydag/results.pkl
