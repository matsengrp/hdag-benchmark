#!/bin/bash

# Uses ...

echo ""
echo "=> Aggregating results..."

set -eu
eval "$(conda shell.bash hook)"
conda activate hdag-benchmark
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

cd ~/hdag-benchmark/data

ignoredir=_ignore
mkdir -p $ignoredir
for clade in $(cat clades.txt); do
    echo $clade
    cladedir=$ignoredir/$clade
    results=$cladedir/results.pkl
    outdir=$cladedir/results
    mkdir -p $outdir

    python ~/hdag-benchmark/support_pipeline_scripts/cli.py agg -i $results -o $outdir
done
