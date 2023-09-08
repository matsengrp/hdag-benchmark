# TODO:
- Describe the purpose of this directory
- Give instructions for how to run experiments


# Description

This folder contains scripts for better understanding the diffused sampling.
- ...

## Approach
We simulate data with varying amounts of hypermutation, run MrBayes on it, and build a coverage plot to evaluate these effects.


# Instructions

1. Run `run_simulations_cluster.sh` to generate a bunch of simulated data.

2. Run `run_inference_cluster.sh <method>` to use the given `method` for support estimation. 

3. To create boxplots for understanding the relationship between parsimony diversity and hypermutation run `python support_pipeline/plotting.py plot_hmut -d /fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation`. The results will be stored under `hypermutation/figures`

4. To create coverage plots for MrBayes, run `bash support_pipeline/scripts/hypermutation_experiment/plot_results.sh mrbayes`