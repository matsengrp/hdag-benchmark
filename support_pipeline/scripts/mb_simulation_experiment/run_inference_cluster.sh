#!/bin/bash

num_trials=5
num_cores=1

method=$1


echo ""
echo "=> Starting inference for "$method"..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

currdir=$PWD
datadir=data/sim_models
cd $datadir
script_dir=$currdir"/support_pipeline/scripts/mb_simulation_experiment"
for clade in $(cat $script_dir/clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade

    # Randomly resolve
    for trial in $(seq $num_trials); do
        # let "trial=$trial_+5" 

        cat $script_dir/sim_params.txt | while read line
        do
            [[ $line =~ ^#.* ]] && continue

            IFS=':' read -a arr <<< $line
            echo "Simulation name: ${arr[0]}"

            resultsdir=$clade/${arr[0]}/$trial/"results/$method"
            mkdir -p $resultsdir

            logfile=$resultsdir/inference_${method}.log

            sbatch -c $num_cores -t 10-0 -J "${arr[0]}|$clade|inference" -o $logfile -e $logfile.err \
            $currdir/support_pipeline/scripts/mb_simulation_experiment/infer_trial_${method}.sh $currdir $clade $trial "${arr[0]}"
        
        done
    done
done