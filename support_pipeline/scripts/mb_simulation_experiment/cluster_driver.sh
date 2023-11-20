#!/bin/bash

num_trials=25
num_cores=1


echo ""
echo "=> Gathering MrBayes parsimony posterior..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark

datadir=$PWD/data/sim_models
script_dir=$PWD"/support_pipeline/scripts/mb_simulation_experiment"

for clade in $(cat $script_dir/clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade

    for trial in $(seq $num_trials); do
        out_dir=$PWD/data/sim_models/$clade/$trial
        mkdir -p $out_dir
        logfile=$out_dir/mb_parsimony_distribution.log

        sbatch -c $num_cores -J "$trial|$clade" -o $logfile -e $logfile.err \
        $script_dir/mb_parsimony_distribution.sh $PWD $clade $trial gamma_10_hmut_50
        
    done
done