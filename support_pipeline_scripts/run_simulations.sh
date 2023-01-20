#!/bin/bash

# Simulates SARS-CoV-2 data with the following procedure
#   1. Subset the big UShER generated MAT for a tree topology
#   2. Randomly resolve all multifurcations
#   3. Use PhastSim to simulate internal node sequences all the way out to leaf nodes
#   4. Store leaf sequences, ancestral sequence, and simulated tree

# Paramaters that determine how many trees to simulate for each clade
num_res=2
num_sim=2
# TODO: makes these parameters for inference and simulation scripts

echo ""
echo "=> Simulating trees..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

mkdir -p ~/hdag-benchmark/data
cd ~/hdag-benchmark/data

# Get sarscov usher tree and reference sequence
bigtree=public-2022-10-01.all.masked.pb.gz
[ -f $bigtree ] || {
    echo "Downloading data..."
    wget -q http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2022/10/01/$bigtree
}
wget -q -O refseq.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MN908947.3&rettype=fasta"
sed -i "1s/.*/>ancestral/" refseq.fasta

# remove newlines from sequence body for convenience
head -1 refseq.fasta > _refseq.fasta
tail -n +2 refseq.fasta | tr -d '\n' >> _refseq.fasta
echo >> _refseq.fasta
mv _refseq.fasta refseq.fasta

for clade in $(cat ../clades.txt); do
    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    # TODO: Once working, send a cluster job off that does all this stuff?
    echo $clade
    cladedir=$clade
    mkdir -p $cladedir
    tree=$cladedir/tree.nwk
    matUtils extract -c $clade -i $bigtree -t $tree # Gets the tree topology for subtree of this clade

    # creates numerified tree newick with name $ntree
    ntreebase=$cladedir/tree
    hdb numerify-taxon-names $tree $ntreebase
    ntree=${ntreebase}.n.nwk

    # Randomly resolve
    for rseed in $(seq $num_res); do
        echo "   $rseed"
        rdir=$clade/"resolved"/$rseed
        mkdir -p $rdir
        rtree=$rdir/rtree.nwk
        
        hdb resolve-multifurcations -i $ntree -o $rtree --resolve-seed $rseed

        for sim in $(seq $num_sim); do
            echo "      $sim"
            let "trial = ($rseed-1) * $num_res + $sim"
            # Converts pair of rseed trial into a single 1-indexed number

            simdir=$clade/$trial/"simulation"
            mkdir -p $simdir

            # Simulate mutations on the tree until we get one without convergent evolution
            for seed in $(seq 10); do
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
                ctreefasta_with_refseq=$simdir/ctree_with_refseq.fasta

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
        done
    done
done
