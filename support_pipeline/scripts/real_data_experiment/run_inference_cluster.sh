#!/bin/bash

num_cores=4
method=$1

echo ""
echo "=> Starting inference for $method..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python


currdir=$PWD
datadir=data/real_data
cd $datadir
script_dir=$currdir"/support_pipeline/scripts/real_data_experiment"
for clade in $(cat $script_dir/clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade
    cladedir=$datadir/$clade

    # for trial in $(seq $num_trials); do
    logfile=$clade"/results/"$method-output-rerun.log
    echo $logfile

    mkdir -p $clade"/results"
    
    sbatch -t 7-0 -c $num_cores -J "$clade|$method|real-data-mb" -o $logfile -e $logfile.err \
    $currdir/support_pipeline/scripts/real_data_experiment/run_$method.sh $currdir $clade
done

