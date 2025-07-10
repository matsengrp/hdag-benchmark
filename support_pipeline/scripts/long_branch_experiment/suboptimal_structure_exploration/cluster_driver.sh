#!/bin/bash

num_trials=25
num_cores=1


echo ""
echo "=> Starting PCM anlaysis..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark

datadir=$PWD/data/sim_models
script_dir=$PWD"/support_pipeline/scripts/suboptimal_structure_exploration"

for clade in $(cat $script_dir/clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade

    for trial in $(seq $num_trials); do
        out_dir=$PWD/data/sub_struct/$clade/$trial
        mkdir -p $out_dir
        logfile=$out_dir/compare_trial.log

        sbatch -c $num_cores -J "$trial|$clade" -o $logfile -e $logfile.err \
        $script_dir/compare_trial.sh $PWD $clade $trial
        
    done
done