#!/bin/bash

# Reconstructs the true sequences using matUtils MAT datastructure

num_cores=1

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

currdir=$PWD

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

cd real_data

ls ../..

echo "In directory $PWD"

script_dir=$currdir"/support_pipeline/scripts/real_data_experiment"
for clade in $(cat $script_dir/clades.txt); do
    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade
    cladedir=$clade
    mkdir -p $cladedir
    tree=$cladedir/tree.nwk
    # Gets the tree topology for subtree of this clade with num mutations on edges
    matUtils extract -c $clade -i ../../$bigtree -t $tree

    # creates numerified tree newick with name $ntree
    ntreebase=$cladedir/tree
    hdb numerify-taxon-names $tree $ntreebase
    ntree=${ntreebase}.n.nwk
    export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
    pb_tree=$cladedir/clade_tree.pb.gz
    matUtils extract -c $clade -i ../../$bigtree -o $pb_tree

    conda activate bte
    python /fh/fast/matsen_e/whowards/hdag-benchmark/support_pipeline/scripts/real_data_experiment/reconstruct_fasta.py $PWD/$cladedir
    conda activate hdag-benchmark
done