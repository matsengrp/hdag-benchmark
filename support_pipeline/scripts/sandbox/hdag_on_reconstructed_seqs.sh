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

datadir=$currdir
trialdir=$datadir/$clade/real/$trial
resultsdir=$trialdir/results
mkdir -p $resultsdir

ctree=$trialdir/collapsed_hdag_tree.nwk
# ctreefasta_with_refseq=$ctree.fasta
# ctreevcf=${ctreefasta_with_refseq}.vcf
# faToVcf $ctreefasta_with_refseq $ctreevcf

dagdir=$resultsdir/historydag
mkdir -p $dagdir
# starttree=$dagdir/starttree.nwk
# echo $dagdir
# echo "()ancestral;" > "$starttree"
# refseqfile=$dagdir/refseq.txt
# head -2 $ctreefasta_with_refseq | tail -n +2 > $refseqfile
# seedtree=$dagdir/seedtree.pb

mkdir -p $dagdir/opt_info
optdag_final=final_opt_dag

# Copying olddag to new dir for now
cp $datadir/$clade/real/results/historydag/final_opt_dag.pb $datadir/$clade/real/$trial/results/historydag/final_opt_dag.pb
cp $datadir/$clade/real/collapsed_hdag_tree.nwk.fasta $datadir/$clade/real/$trial/collapsed_hdag_tree.nwk.fasta


# TODO: Uncomment this when you want to search for new trees ---v
# $ctreevcf contains the ancestral sequence and all the other simulated sequences
# echo "===> create tree with UShER..."
# usher-sampled -v $ctreevcf -t $starttree -o $seedtree  --optimization_minutes=0 -d $dagdir/opt_info

# echo "===> lusher optimizing..."
# log_prefix=$dagdir/opt_info/optimization_log
# python ../support_pipeline/inference.py larch_usher -i $seedtree -r $refseqfile -c 4000 -o $dagdir -l $log_prefix -f $optdag_final

# python /fh/fast/matsen_e/whowards/hdag-benchmark/support_pipeline/scripts/sandbox/save_unif_hdag_tree.py $trial
optdag_final='/fh/fast/matsen_e/whowards/hdag-benchmark/data_new_sim/AY.108/real/results/historydag/complete_opt_dag.pb'
python ../support_pipeline/inference.py save_supports -m "hdag" -t $ctree -i $optdag_final -o $dagdir/results.pkl


echo ""
echo ""
# sbatch -c 4 -J "real inf" \
# -o /fh/fast/matsen_e/whowards/hdag-benchmark/data_new_sim/AY.108/real/results/inference_historydag.log \
# -e /fh/fast/matsen_e/whowards/hdag-benchmark/data_new_sim/AY.108/real/results/inference_historydag.err \
# /fh/fast/matsen_e/whowards/hdag-benchmark/support_pipeline/scripts/sandbox/hdag_on_reconstructed_seqs.sh $PWD AY.108 2

