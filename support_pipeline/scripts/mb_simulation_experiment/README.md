# Description

We are curious about why MrBayes is having a hard time producing good coverage plots.
Our hypothesis is that this is caused by a model mismatch.
Our simulations use the UNREST substitution model (a non-reversible model) and hypermutation, whereas MrBayes uses GTR and doesn't consider hypermutation.

# Approach

To better understand the effects of model mismatch on MrBayes coverage, we ran a bunch of different simulations using increasingly complex models, from JC with no rate variation to UNREST with hypermutation.

# Instructions

1. Run `bash support_pipeline/scripts/mb_simulation_experiment/run_simulations_cluster.sh` to generate simulated data that is increasingly realistic (varies in substitituion model, site invaraince, rate variation, amount of hypermutation). This script uses `sim_params.txt` to determine the phastSim parameters (e.g., substitution model, hypermutation, gamma distributed rate variation, etc.)

2. Run `run_inference_cluster.sh <method>` to use the given `method` for support estimation. 





