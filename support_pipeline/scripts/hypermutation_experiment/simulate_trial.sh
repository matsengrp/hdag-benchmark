#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"

conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

# This should always be the full path to data directory
curr_dir=$1
rseed=$2
simdir=$3
rtree=$4
hypermut_prob=$5
hypermut_rate=$6

echo $1
echo $2
echo $3
echo $4
echo $5
echo $6

# NOTE: Should always be data directory
cd $curr_dir

# Simulate mutations on the tree until we get one without convergent evolution
for seed in $(seq 100); do
    # Seed is also printed to sim.log
    let "s = $seed + ($rseed-1) * 100"
    
    phastSim --outpath $simdir/ --seed $s --createFasta --createInfo \
            --createNewick --createPhylip --treeFile $rtree \
            --scale 0.00003344 \
            --categoryProbs 0.25 0.25 0.25 0.25 --categoryRates 0.1 0.5 1.0 2.0 \
            --reference refseq.fasta --eteFormat 1 \
            --hyperMutProbs $hypermut_prob --hyperMutRates $hypermut_rate
            # --mutationRate JC69
    
    simtree=$simdir/sars-cov-2_simulation_output.tree
    simfasta=$simdir/sars-cov-2_simulation_output.fasta

    ctree=$simdir/collapsed_simulated_tree.nwk
    ctreefasta=${ctree}.fasta
    ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta  # Stores leaf sequences
    
    # produces ctreefasta_with_refseq which doesn't contain reference sequence
    # throws error message if tree exhibits convergent evolution
    if hdb collapse-tree $simtree $simfasta $ctree; then

        cp refseq.fasta $ctreefasta_with_refseq
        cat $ctreefasta >> $ctreefasta_with_refseq
        # This will treat the first sequence (`ancestral`) in the vcf as the
        # reference, but not include it as a sample in the vcf.
        ctreevcf=${ctreefasta_with_refseq}.vcf
        faToVcf $ctreefasta_with_refseq $ctreevcf

        break
    fi
done

# TODO: Should we be looking at the resolved tree???
hdb get-tree-stats -s $simdir

# Uncollapse tree
re_resolved_tree=$simdir/resolved_output.nwk
hdb resolve-multifurcations -i $ctree -o $re_resolved_tree --resolve-seed 1 --branch-len-model num-muts
