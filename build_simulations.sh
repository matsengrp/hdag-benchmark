#!/bin/bash
set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

cd data
# Get sarscov usher tree and reference sequence
bigtree=public-2022-10-01.all.masked.pb.gz

[ -f $bigtree ] || {
    echo "Downloading data..."
    wget -q http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2022/10/01/$bigtree
}
wget -q -O refseq.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MN908947.3&rettype=fasta"
sed -i "1s/.*/>ancestral/" refseq.fasta

# ## # For choosing clades of appropriate sizes:
# ## matUtils summarize -A -i $bigtree
# ## python ../choose_clades.py > clades.txt
ignoredir=_ignore
mkdir -p $ignoredir
for clade in $(cat clades.txt); do
    # send a cluster job off that does all this stuff?

    cladedir=$ignoredir/$clade/
    mkdir -p $cladedir
    tree=$cladedir/tree.nwk
    matUtils extract -c $clade -i $bigtree -t $tree

    # creates numerified tree newick with name $ntree
    ntreebase=$cladedir/tree
    hdb-numerify-taxon-names.sh $tree $ntreebase
    ntree=${ntreebase}.n.nwk

    # creates output files $cladedir/sars-cov-2_simulation_output.tree an
    # extended newick containing mutations and
    # $cladedir/sars-cov-2_simulation_output.fasta containing leaf sequences
    # (among other files):
    phastSim --outpath $cladedir --seed 9 --createFasta --createInfo \
             --createNewick --createPhylip --treeFile $ntree \
             --scale 0.00004 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 \
             --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon \
             --reference refseq.fasta --eteFormat 1
    simtree=$cladedir/sars-cov-2_simulation_output.tree
    simfasta=$cladedir/sars-cov-2_simulation_output.fasta

    ctree=$cladedir/collapsed_simulated_tree.nwk
    ctreefasta=${ctree}.fasta
    ctreefasta_with_refseq=$cladedir/ctree_with_refseq.fasta

    # produces ctreefasta_with_refseq which contains reference sequence
    hdb collapse-tree $simtree $simfasta $ctree
    cp refseq.fasta $ctreefasta_with_refseq
    cat $ctreefasta >> $ctreefasta_with_refseq



    ##### Start inference: build dag with dnapars
    dagdir=$cladedir/historydag
    echo $dagdir
    rm -rf $dagdir
    find-dnapars-trees.sh -f $ctreefasta_with_refseq -n 1 -o $dagdir -r ancestral
done
