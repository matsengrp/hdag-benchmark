#!/bin/bash

source "/n/fs/ragr-research/users/wh8114/anaconda3/bin/activate"

num_trials=5
num_cores=5

method=historydag

# NOTE: Need to be in larch-env when running this

echo ""
echo "=> Starting inference for "$method"..."

# set -eu
# eval "$(conda shell.bash hook)"
conda activate hdag-benchmark

# Assumes that we are in the hdag-benchmark directory
currdir=$PWD
datadir=data/sim_models
cd $datadir
script_dir=$currdir"/support_pipeline/scripts/long_branch_experiment"
for clade in $(cat $script_dir/clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade

    for trial in $(seq $num_trials); do
        # for branch_multiplier in 2 4 8 16 32 64 128 256; do
        for branch_multiplier in 1 16 256; do

            trial_=$((trial))

            trialdir=$currdir/$datadir/$clade/$trial_/branch_multiplier_$branch_multiplier
            resultsdir=$trialdir/$method

            echo "Making $resultsdir"
            mkdir -p $resultsdir

            logfile=$resultsdir/inference_${method}.log

            sbatch -c $num_cores -t 6-0 -J "$trial_|$branch_multiplier" -o $logfile -e $logfile.err \
            $currdir/support_pipeline/scripts/long_branch_experiment/infer_trial_${method}.sh $currdir $trialdir

        done
    done
done