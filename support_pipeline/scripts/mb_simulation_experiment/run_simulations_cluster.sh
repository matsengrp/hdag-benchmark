#!/bin/bash


# Paramaters that determine how many trees to simulate for each clade
num_res=25
num_sim=1
num_cores=1

echo ""
echo "=> Simulating trees..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

# Assume we are in one directory above `data`
currdir=$PWD

mkdir -p data/sim_models
cd data/sim_models

# Get sarscov usher tree and reference sequence
bigtree=../public-2022-10-01.all.masked.pb.gz
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
script_dir=$currdir"/support_pipeline/scripts/mb_simulation_experiment"
for clade in $(cat $script_dir/clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade
    cladedir=$clade
    mkdir -p $cladedir
    tree=$cladedir/tree.nwk
    # Gets the tree topology for subtree of this clade with num mutations on edges
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

        cat $script_dir/sim_params.txt | while read line
        do
            [[ $line =~ ^#.* ]] && continue

            IFS=':' read -a arr <<< $line
            echo "Sim: ${arr[0]} Params: ${arr[1]}";
            
            simdir=$clade/${arr[0]}/$rseed/"simulation"
            mkdir -p $simdir

            sbatch -c $num_cores -J "${arr[0]}|$clade|sim" -o $simdir/sim.log \
            $currdir/support_pipeline/scripts/mb_simulation_experiment/simulate_trial.sh $PWD $rseed $simdir $rtree "${arr[1]}"
        
        done
    done
done