#!/bin/bash

num_trials=25
num_cores=1


echo ""
echo "=> Starting scaling anlaysis..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark

datadir=$PWD/data/sim_models
script_dir=$PWD"/support_pipeline/scripts/rf_dist_scaling_experiment"

for clade in $(cat $script_dir/clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade

    for trial in $(seq $num_trials); do
        out_dir=$PWD/data/rf_scaling/$clade/$trial
        mkdir -p $out_dir
        logfile=$out_dir/rf_scaling

        sbatch -c $num_cores -J "$trial|$clade" -o $logfile"_min_trim.log" -e $logfile"_min_trim.log.err" \
        $script_dir/rf_scaling.sh $clade $trial 'min_trim'
        sbatch -c $num_cores -J "$trial|$clade" -o $logfile"_naive.log" -e $logfile"_naive.log.err" \
        $script_dir/rf_scaling.sh $clade $trial 'naive'
        
    done
done