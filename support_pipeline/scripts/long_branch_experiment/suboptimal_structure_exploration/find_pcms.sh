#!/bin/bash

source /n/fs/ragr-research/users/wh8114/anaconda3/bin/activate hdag-benchmark

# This should always be the full path to data directory
clade=$1
trial=$2
branch_multiplier=$3

echo "clade name: $1"
echo "trial:      $2"
echo "branch:     $3"

results_dir="/n/fs/ragr-research/users/wh8114/prev/hdag-benchmark/data/sub_struct/results"
trial_dir=$results_dir/$clade/$trial/branch_multiplier_$branch_multiplier
seedtree=$trial_dir/sim_history.pb

echo
echo "==> saving simulated tree as history..."
python support_pipeline/scripts/long_branch_experiment/suboptimal_structure_exploration/save_sim_history.py $clade $trial $branch_multiplier

source /n/fs/ragr-research/users/wh8114/anaconda3/bin/activate larch-env

echo
echo "==> lusher optimizing..."
log_prefix=$trial_dir/opt_info/optimization_log_sim
mkdir -p $log_prefix

# TODO: Trying no --trim
# NOTE: This doesn't work well....

larch-usher \
-i $seedtree \
-o $trial_dir/sim_dag.pb \
--move-coeff-nodes 1 \
--move-coeff-pscore 10 \
--sample-any-tree \
-l $log_prefix \
-c 5 \
--input-format dag-pb \
--output-format dag-pb

echo
echo "==> finding PCMs..."
source /n/fs/ragr-research/users/wh8114/anaconda3/bin/activate hdag-benchmark
python support_pipeline/scripts/long_branch_experiment/suboptimal_structure_exploration/find_pcms.py $clade $trial $branch_multiplier
