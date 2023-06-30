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
datadir=data
cd $datadir

for clade in $(cat ../clades.txt); do

    # if value of $var starts with #, ignore it
    [[ $clade =~ ^#.* ]] && continue

    echo $clade
    cladedir=$datadir/$clade"_"

    # for trial in $(seq $num_trials); do
    logfile=$clade"_/results/"$method-output.log
    echo $logfile

    mkdir -p $clade"_/results"
    
    sbatch -t 7-0 -c $num_cores -J "$clade|$method|real-data-mb" -o $logfile -e $logfile.err \
    $currdir/support_pipeline/scripts/real_data_experiment/run_$method.sh $currdir $clade

done

