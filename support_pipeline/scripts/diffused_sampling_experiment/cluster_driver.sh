#!/bin/bash


num_trials=25
num_cores=1

# NOTE: Don't include .sh in the script name. Script must be an executable that takes
#       the current working directory, clade name, and trial number
script_name=$1

echo ""
echo "=> Running $script_name with $num_trials trials on $num_cores cores..."

set -eu
eval "$(conda shell.bash hook)"
conda activate base

datadir=$PWD/data/sim_models
script_dir=$PWD"/support_pipeline/scripts/diffused_sampling_experiment"

for clade in $(cat $script_dir/clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade

    for trial in $(seq $num_trials); do
        out_dir=$PWD/data/sim_models/$clade/gamma_10_hmut_50/$trial/results/diff_historydag
        # NOTE: out_dir should already exist. If not, it's probably incorrectly formatted
        # mkdir -p $out_dir
        logfile=$out_dir/${script_name}.log

        sbatch -c $num_cores -J "$trial|$clade" -o $logfile -e $logfile.err \
        $script_dir/cluster_script.sh $script_name $clade $trial
        
    done
done