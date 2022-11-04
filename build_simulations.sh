#!/bin/bash
set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
larch_usher_exec=~/larch-usher3/larch/build/larch-usher

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

# ## # For choosing clades of appropriate sizes:
# ## matUtils summary -A -i $bigtree
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

    # Simulate mutations on the tree until we get one without convergent
    # mutations
    for seed in $(seq 10); do
        phastSim --outpath $cladedir --seed $seed --createFasta --createInfo \
                 --createNewick --createPhylip --treeFile $ntree \
                 --scale 0.00004 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 \
                 --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon \
                 --reference refseq.fasta --eteFormat 1
        simtree=$cladedir/sars-cov-2_simulation_output.tree
        simfasta=$cladedir/sars-cov-2_simulation_output.fasta

        ctree=$cladedir/collapsed_simulated_tree.nwk
        ctreefasta=${ctree}.fasta
        ctreefasta_with_refseq=$cladedir/ctree_with_refseq.fasta

        # produces ctreefasta_with_refseq which doesn't contain reference sequence
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



    ##### Start inference: build dag with larch-usher
    dagdir=$cladedir/historydag
    mkdir -p $dagdir
    starttree=$dagdir/starttree.nwk
    echo "();" > $starttree
    refseqfile=$dagdir/refseq.txt
    head -2 $ctreefasta_with_refseq | tail -n +2 > $refseqfile
    seedtree=$dagdir/seedtree.pb
    
    usher-sampled -v $ctreevcf -t $refseqfile -o $seedtree --optimization_minutes=0
    $optimized_dag=$dagdir/opt_dag.pb
    $larch_usher_exec -i $seedtree -r $refseqfile -c 100 -o $optimized_dag
done
