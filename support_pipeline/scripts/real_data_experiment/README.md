# Description

This folder contains scripts for running a particular type of experiment.
Rather than simulate data, we use the real data (as reconstructed from MAT datastructure and ancestral sequence), and run MrBayes on the resulting sequences.
We then sample a tree from the MB posterior, treat that as the true tree, and use CA plots to evaluate our method wrt to this tree.

We are shifting to this evaluation procedure because it was hard to simulate realistic data.
Specifically, in every simulation, there was little MP diversity (number of MP trees and types of nodes) on that data.

Included scripts are variations on simulation, inference, and plotting scripts from parent directory. 

# Instructions

1. Run `reconstruct_fasta.sh` to create the leaf sequences

2. Run MrBayes inference with `cluster_driver.sh mrbayes`
    - Sample topology from MrBayes tree file (sample line after burnin)

3. Run HistoryDAG inference with `cluster_driver.sh historydag`

- At this point, `data/AY.108_/results/<method>/results.pkl` will contain results lists, where `<method>` is historydag or mrabyes

4. Plot CA results with `python support_pipeline/plotting.py aggregate_results`
    