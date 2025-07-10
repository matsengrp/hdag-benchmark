#!/bin/bash

source "/n/fs/ragr-research/users/wh8114/anaconda3/bin/activate"
conda activate larch

num_cores=5
intermediate_paths="/n/fs/ragr-research/users/wh8114/prev/hdag-benchmark/data/sim_models/AY.34.2/intermediate_file_paths.txt"

for dag_path in $(cat $intermediate_paths); do

    # if value of $var starts with #, ignore it
    [[ $dag_path =~ ^#.* ]] && continue

    results_dir="$dag_path"_dir
    echo "Making $results_dir"
    mkdir -p $results_dir
    logfile=$results_dir/trimming.log
    sbatch -c $num_cores -t 6-0 -o $logfile -e $logfile.err \
    /n/fs/ragr-research/users/wh8114/prev/hdag-benchmark/support_pipeline/scripts/long_branch_experiment/trim_dag.sh \
    $dag_path $results_dir

    # Use this to copy into main results dir for PCM finding
    # cp $dag_path"_dir/final_opt_dag_trimmed.pb" $dag_path"_dir/.."

done