# hdag-benchmark

This repository contains scripts for running experiments and generating figures for an in-review paper *Densely sampled phylogenies frequently deviate from maximum parsimony in simple and local ways* by William Howard-Snyder, Will Dumm, Mary Barker, Ognian Milanov, Claris Winston, David H. Rich, and Frederick A Matsen IV.

## Dependencies

Clone this repository, and `cd` into it.
Install dependencies via:

    conda env create -f environment.yml
    conda activate hdag-benchmark
    pip install -e .

We also make heavy use of the [historyDAG repository](https://github.com/matsengrp/historydag).
The branch that we're using currently is `whs-diffused-madag-sampler`.
You can access it by cloning the repository and running `pip install -e .` once inside.
We also use the C++ implementation of the history sDAG, Larch, from [this repo](https://github.com/matsengrp/larch).
You can access it by cloning the repository and running the build instructions.

Many of the scripts described below use assume access to a cluster and attempt to access it via the slurm commands.

## Main Results Replication: _Densely sampled phylogenies frequently deviate from maximum parsimony in simple and local ways_

Here we describe how to use this codebase to reproduce our main results.
All scripts should be run from immediately within the `hdag-benchmark` directory.

## Simulations
We simulate data by running

    bash support_pipeline/scripts/mb_simulation_experiemnt/run_simulations_cluster.sh

This requires that each line of `sim_params.txt`corresponds to the rate variation commands we provide for phastSim.
The command that we used to simulate data for our main results was

    gamma_10_hmut_50:--alpha 0.1 --hyperMutProbs 0.01 --hyperMutRates 50.0

which populates `data/sim_models/<clade>/<trial>/gamma_10_hmut_50/simulation` with the simulated data for each clade in `clades.txt` and trial.

## MP Tree Search
We produce large sets of histories stored as Mutation Annotated historyDAGs by running

    bash support_pipeline/scripts/mb_simulation_experiemnt/run_inference_cluster.sh

The resulting protobufs will be stored at this filepath for each clade and trial

    data/sim_models/<clade>/<trial>/gamma_10_hmut_50/results/historydag/final_opt_dag_trimmed.pb

and used for replicating the following results.

## Proportion due to PCM (Figure 4)
To generate the number of differences and proportion that are due to PCMs run

    bash support_pipeline/scripts/suboptimal_cluster_exploration/cluster_driver.sh

which runs the `compare_trial.py` script for each clade and trial. To plot the swarmplot from Figure 4 of our paper, run

    bash support_pipeline/plotting.py plot_sub_stats -d ~/hdag-benchmark/data/sub_struct

This will produce and save the swarmplot and histogram at `hdag-benchmark/data/sub_struct/figures/proportion_duplicate_swarm_per_clade.pdf`.


## RF-Scaling (Figure 3a)

To generate the time data for finding the minimum RF-distance tree run

    bash support_pipeline/scripts/rf_dist_scaling_experiment/cluster_driver.sh

which runs the `rf_scaling.py` script for each clade and trial. To plot the scatterplot from Figure 4 of our paper, adjust the `datapath` variable in the `plot_rf_scaling()` method of `support_pipeline/plotting.py` to point to the `data/rf_scaling` directory and run

    bash support_pipeline/plotting.py plot_rf_scaling

This will save the scatterplot to `hdag-benchmark/data/rf_scaling/figures/scaling.pdf`.



## Directory Description

The rest of this README describes the rest of the code that is not used in our most recent paper.

### Scripts

The `support_pipeline` directory should have the following structure:

    ├── data
    ├── hdb
    └── support_pipeline
        ├── simulation.py
        ├── inference.py
        ├── plotting.py
        ├── plot_utils.py
        ├── utils.py
        └── scripts
            ├── diffused_sampling_experiment
            ├── hypermutation_experiment
            ├── mb_simulation_experiment
            ├── suboptimal_structure_exploration
            ├── real_data_experiment
            ├── plot_results.sh
            ├── python_replace.py
            ├── run_inference_cluster.sh
            ├── run_simulations_cluster.sh
            ├── infer_trial_historydag.sh
            ├── infer_trial_mrbayes.sh
            ├── infer_trial_random.sh
            ├── infer_trial_ufboot.sh
            └── simulate_trial.sh

The files `simulation.py`, `inference.py`, and `plotting.py` contain cli scripts for simulating data, inferring support/posteriors, and plotting results, respectively.
These python cli scripts are bundled up into bash files that can be sent to the cluster by cluster drivers (bash scripts that have `cluster` in their name).
For example, running

    bash support_pipeline/scripts/run_inference_cluster.sh historydag

executes `infer_trial_historydag.sh` for each clade and trial you've specified, and sends these instances to the cluster.

The directories under scripts that end in experiment (or exploration) contain scripts that we used to investigate a particular phenomenon (e.g., the relationship between hypermutation and parsimony diversity) that isn't necessarily directly related to support inference.
Often, these will have there own cluster drivers, and simulation/inference scripts that are minor variations of those in the `scripts` directory.


### Data + Results

Once you run a simulation script, it will populate the `data` directory.

Currently, we are simulating using hypermutation and gamma rate variation. You can simulate this data by running

    bash support_pipeline/scripts/mb_simulation_experiemnt/run_simulations_cluster.sh

Note that the cluster driver is within the `mb_simulation_experiment`, NOT the one directly under the scripts directory.
Running this command will produce the following directory structure:

    ├── data
    |   └── sim_models
    |       └── <clade-name>
    |           └── gamma_10_hmut_50
    |               ├── <trial-number>
    |               |   ├── figures
    |               |   ├── results
    |               |   └── simulation
    |               └── figures
    |
    ├── hdb
    └── support_pipeline

The simulated data for the given clade-trial will be stored in the `data/sim_models/<clade-name>/gamma_10_hmut_50/<trial-number>/simulation` directory, and figures/results for that trial will be stored in the figure/results directory, respectively.

The `figures` directory outside the trial directory is used to store figures that are aggregated over all the trials for that clade (e.g., coverage analysis plot for all simulations that use the A.2.5 UShER-clade for their topology).

## Experiments

Most of the results/data for these experiemnts are not pushed to the repository.
You can find them by going to my fast directory `fh/fast/whowards/hdag-benchmark` and appending the given file path stated below.
You can read more about each experiment below or by clicking their links:

- [`hypermutation_experiment`](https://github.com/matsengrp/hdag-benchmark/tree/main/support_pipeline/scripts/hypermutation_experiment):
Used to observe how hypermutation impacts a method's coverage performance. Results + data stored in `data/hypermutation`.

- [`real_data_experiment`](https://github.com/matsengrp/hdag-benchmark/tree/main/support_pipeline/scripts/real_data_experiment):
Evaluates methods of support estimation by running MrBayes on real data, sampling a tree from the MrBayes posterior, and treating it as the true tree when producing coverage plot. Results + data stored in `data/real_data` and figures at `data/figures`.

- [`mb_simulation_experiment`](https://github.com/matsengrp/hdag-benchmark/tree/main/support_pipeline/scripts/mb_simulation_experiment):
Evaluates MrBayes support/parsimony posterior coverage at varying levels of model mismatch.
Also, used these scripts investigate Parsimony Diversity.
Results + data stored at `data/sim_models`.

- [suboptimal_structure_exploration](https://github.com/matsengrp/hdag-benchmark/tree/main/support_pipeline/scripts/suboptimal_structure_exploration):
Explores what sorts of unparsimonious structures are present in the true trees.
Results stored at `data/sub_struct`.

- [diffused_sampling_experiment](https://github.com/matsengrp/hdag-benchmark/tree/main/support_pipeline/scripts/diffused_sampling_experiment):
Used to test and visualize the diffused sampler.
Also, contains scripts for identifying when a mutation in the MP tree is likely to be a PCM. 
Results stored in the `diff_historydag` directory of the results section of `sim_models` (e.g., `data/sim_models/A.2.5/gamma_10_hmut_50/figures/diff_historydag_mut_rates`).

The first two are somewhat irrelevant now. So they can be safely ignored.
The experiments in `mb_simulation_experiment` are not so important, but we currently use the simulation scripts to simulate data that has hypermutation and gamma rate variation.
This simulated data is used by the last two experiments, which are more relevant to what we're working on now.

All bash scripts should be run from the hdag-benchmark directory, unless otherwise specified.


## Simulation

The general pipeline for simulating data you can run a batch of simulations with 

    bash <path-to-simulation-script>/run_simulations_cluster.sh

where `<path-to-simulation-script>` depends on whether you are running the data simulations for an experiment or not (e.g., `support_pipeline/scripts/mb_simulation_experiment`). Regardless, this method will extract tree topologies from the big UShER MAT, and simulates sequences for each clade using PhastSim on those topologies.

If using the cluster driver in the `scripts` directory, these simulations can be found at `data/<clade-name>/<trial-number>/simulation`.
Otherwise, they will be stored at a `simulation` directory that has a similar, but slightly different path that may depend on the experiment it came from.
In this directory, `ctree_with_refseq.fasta` contains the leaves of the simulated tree and the ancestral sequence and `collapsed_simulated_tree.nwk` contains the simulated tree.

## Inference

After running simulations, you can infer supports using 

    bash <path-to-inference-script>/run_inference_cluster.sh <method-name>

This will produce a `results.pkl` in the directory

    data/<clade-name>/<trial-number>/results/<method>

for each clade and trial. This is a pickled list of tuples `(clade, estimate-support, in-tree)`, where

1. `clade` is the set of all taxa under a particular node.
2. `estimate-support` is the inferred support for that node.
3. `in-tree` is a boolean that indicates whether that node is in the true tree.

This directory also includes other information about the inference such as a log of trees for MrBayes and the optimized dag.
