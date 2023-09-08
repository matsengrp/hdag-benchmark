# Description

This folder contains scripts for exploring what sorts of unparsimonious structures are present in the simulated trees and what sorts of mutations give rise to these structures. The hope is that these will inform a correction of branch support.


# Instructions

*Assumes `data/sim_models` has been populated with 25 trials of `gamma_10_hmut_50` simulated data*

1. Run `bash support_pipeline/scripts/suboptimal_structure_exploration/cluster_driver.sh`.
This step generates PDFs in `data/sub_struct/<clade>/<trial>` of the simulated tree and the most similar MP tree and identifies where structural differences occur.
It also produces a log file `compare_trial.log` for each trial, that states the number of differing nodes and the number that are due to a PCM.

2. Run `python support_pipeline/plotting.py plot_sub_stats -d /fh/fast/matsen_e/whowards/hdag-benchmark/data/sub_struct`.
This step parses the information from the log file and presents it as a swarm plot that gets stored at `sub_struct/figures`.
