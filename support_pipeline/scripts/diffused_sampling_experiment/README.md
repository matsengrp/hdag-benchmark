# Description

This folder contains scripts for and plotting parsimony posterior that it induces.


# Instructions

*Assumes `data/sim_models` has been populated with 25 trials of `gamma_10_hmut_50` simulated data*

Run `bash support_pipeline/scripts/diffused_sampling_experiment/cluster_driver.sh <name-of-script>` to run the script with the given name (in this directory) for each clade and trial on the cluster.

The two options for `<name-of-script>` are currently,

1. `diff_hdag_parsimony_posterior.py`: This samples a bunch of trees from the diffused dag distribution, and stores the parsimony scores as a pickled counter.
These can be plotted by running

        python support_pipeline/plotting.py plot_pars_distribution \
        -d /fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models \
        -i diff_historydag_mut_rates/pars_distr.pkl

and will be stored in that clade's `figures` directory (e.g., `data/sim_models/A.2.5/gamma_10_hmut_50/figures/diff_historydag_mut_rates/parsimony_distribution_relative_change.png`)

2. `prob_pcm.py`: For each mutation in the simulated tree, stores a pickled list which contains the number of times it occurred in a PCM, the number of times it occurred, and its mutation rate.
These can be plotted by running

        python support_pipeline/plotting.py examine_pcm_prob \
        -d /fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models 

and will be stored in that clade's `figures` directory