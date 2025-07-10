#!/bin/bash

# Script for running larch on specified seed tree (from UShER) and reference sequence.
# Expects these files to be in the curr_dir directory. Will send all larch output to curr_dir.

# Change these as needed
curr_dir=/n/fs/ragr-research/users/wh8114/prev/hdag-benchmark/temp    # Directory you are working in
larch_dir=/n/fs/ragr-research/users/wh8114/prev/larch                 # Path to larch directory

log_prefix=$curr_dir/opt_info/optimization_log
mkdir -p $log_prefix

seed_tree_path=$curr_dir/seedtree.pb_tree
refseq_path=$curr_dir/refseq.txt


# Replace this with the path to your conda directory
source /n/fs/ragr-research/users/wh8114/anaconda3/bin/activate larch

$larch_dir/build/bin/larch-usher \
-i $seed_tree_path \
-o $curr_dir/final_opt_dag_trimmed.pb \
-r $refseq_path \
--move-coeff-nodes 1 \
--move-coeff-pscore 5 \
--sample-any-tree \
--trim \
-l $log_prefix \
--inter-save 10 \
-c 100 \
--input-format tree-pb \
--output-format dag-pb