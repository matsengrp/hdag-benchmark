#!/bin/bash

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
# NOTE: This file must have executable permissions to be used with simulation pipeline

curr_dir=$1
rseed=$2
sim=$3
simdir=$4
rtree=$5

cd $curr_dir

# Simulate mutations on the tree until we get one without convergent evolution
for seed in $(seq 50); do
    # Seed is also printed to sim.log
    let "s = $seed + ($sim-1) * 10 + ($rseed-1) * 100"
    echo "seed $s"
    phastSim --outpath $simdir/ --seed $s --createFasta --createInfo \
            --createNewick --createPhylip --treeFile $rtree \
            --scale 0.00004 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 \
            --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon \
            --reference refseq.fasta --eteFormat 1
    simtree=$simdir/sars-cov-2_simulation_output.tree
    simfasta=$simdir/sars-cov-2_simulation_output.fasta

    ctree=$simdir/collapsed_simulated_tree.nwk
    ctreefasta=${ctree}.fasta
    ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta  # Stores leaf sequences

    # You made a bunch of new nodes when resolving so your tree.mapping
    # is not completely correct any more
    
    # produces ctreefasta_with_refseq which doesn't contain reference sequence
    # throws error message if tree exhibits convergent evolution
    if hdb collapse-tree $simtree $simfasta $ctree; then
        break
    fi
done
cp refseq.fasta $ctreefasta_with_refseq
cat $ctreefasta >> $ctreefasta_with_refseq
# This will treat the first sequence (`ancestral`) in the vcf as the
# reference, but not include it as a sample in the vcf.
ctreevcf=${ctreefasta_with_refseq}.vcf
faToVcf $ctreefasta_with_refseq $ctreevcf