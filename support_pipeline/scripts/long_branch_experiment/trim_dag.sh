#!/bin/bash

source "/n/fs/ragr-research/users/wh8114/anaconda3/bin/activate"
conda activate larch

dag_path=$1
results_dir=$2

/n/fs/ragr-research/users/wh8114/prev/larch/build/bin/larch-dagutil \
    -i $dag_path \
    -o $results_dir/final_opt_dag_trimmed.pb \
    --trim \
    --input-format dag-pb \
    --output-format dag-pb