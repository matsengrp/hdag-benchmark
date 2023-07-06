#!/bin/bash

# Uses usher (and other methods later) to create a starting tree on the simulated data, then uses larch-usher to optimize
# This script is intended to be run after `run_simulations.sh`
# Should output file containing node supports for every simulated dataset

# Paramaters that determine how many trees to simulate for each clade
num_res=3
num_sim=1
num_trials=$num_res*$num_res
num_cores=2
trial_file="pars_div_trials"
method=$1

let "num_trials = $num_sim * $num_res"

echo ""
echo "=> Starting inference for "$method"..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python


currdir=$PWD
datadir=data/hypermutation
cd $datadir
script_dir=$currdir"/support_pipeline/scripts/hypermutation_experiment"

for clade in $(cat $script_dir/clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade
    cladedir=$datadir/$clade

    for hmut in $(cat $script_dir/hypermutation.txt); do

        for trial in $(seq $num_trials); do

            resultsdir=$clade/$hmut/$trial/"results"
            mkdir -p $resultsdir

            logfile=$resultsdir/inference_$method.log

            sbatch -t 6-0 -c $num_cores -J "$trial|$clade|$hmut" -o $logfile -e $logfile.err \
            $currdir/support_pipeline/scripts/hypermutation_experiment/infer_trial_$method.sh $currdir $clade $trial $hmut
        
        done
    
    done
done

