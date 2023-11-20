# Description

This folder contains scripts for gathering data comparing the time scaling of our MinTrim to a naive minimum distance approach.


# Instructions

*Assumes `data/sim_models` has been populated with 25 trials of `gamma_10_hmut_50` simulated data for each clade*

1. Run `bash support_pipeline/scripts/rf_dist_scaling_experiment/cluster_driver.sh`.
This step saves statistics about the number of histories in the historydag and the time it takes to compute the minimum rf distance tree using naive search and the min_trim historydag method. Each directory of the form `data/rf_scaling/results/<clade>/<trial>` will be populated

2. Run `python support_pipeline/plotting.py plot_rf_scaling`. This step creates a scatter plot summarizing all of the scaling results.
