#!/bin/bash

# Simulates SARS-CoV-2 data with the following procedure
#   1. Subset the big UShER generated MAT for a tree topology
#   2. Randomly resolve all multifurcations
#   2.5 Give added edges relatively small branch lengths
#   3. Use PhastSim to simulate internal node sequences all the way out to leaf nodes
#   4. Collapse tree so that all edges have at least one mutation
#   5. Store leaf sequences, ancestral sequence, and simulated tree

# Paramaters that determine how many trees to simulate for each clade
num_res=1
num_sim=1
num_cores=1

echo ""
echo "=> Simulating trees..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

mkdir -p data
cd data

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

    echo $clade
    cladedir=$clade
    mkdir -p $cladedir
    tree=$cladedir/tree.nwk
    # Gets the tree topology for subtree of this clade with num muations on edges
    matUtils extract -c $clade -i $bigtree -t $tree

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
        
        hdb resolve-multifurcations -i $ntree -o $rtree --resolve-seed $rseed --branch-len-model num-muts

        for sim in $(seq $num_sim); do
            echo "      $sim"
            let "trial = ($rseed-1) * $num_res + $sim"
            # Converts pair of rseed trial into a single 1-indexed number

            simdir=$clade/$trial/"simulation"
            mkdir -p $simdir

            sbatch -c $num_cores -J "$trial|$clade|sim" -o $simdir/sim.log \
            ./../support_pipeline/scripts/simulate_trial.sh $PWD $rseed $sim $simdir $rtree
        done
    done
done
