# hdag-benchmark

## Dependencies

Install dependencies via:

    conda env create -f environment.yml
    conda activate hdag-benchmark


## Install

    pip install -e .


## Additional Requirements

TODO: Add requirements for installing MrBayes v3.2.7a from source (without Beagle or MPI)

We also make use of `beast2-xml`, a python package for loading fasta files into XML format.
Unfortunately, this package is incompatible with `protobuf`.
So, we create a new conda environment to switch to when loading the XML files as follows:

    conda create -n beast-xml python=3.10
    conda activate beast-xml
    pip install beast2-xml

We also make ues of some packages in the most recent version of history dag (specifically, the partial-sankoff branch).
To use these, clone the history dag repository and run
    
    conda activate hdag-benchmark
    pip install -e <path/to/hdag/repo>
    
Make sure you are not in the `beast2-xml` environment, as this could break it.

## Support Experiments

Currently, our scripts for running support experiments is located in `support_pipeline_scripts`
(this might change soon and they'll go in `hbd/scripts`).
These scripts assume that you are running them from the `hdag-benchmark` directory.
Edit the `clades.txt` file to determine which PANGO clades will be simulated/analyzed.
If a `#` is put in front of a clade name (no space), it will be ignored.
The `*_cluster.sh` scripts assume access to a cluster with slurm workload manager.

### Simulation

You can run a batch of simulations with 

    bash support_pipeline_scripts/run_simulations_cluster.sh

This extracts tree topologies from teh big UShER MAT, and simulates sequences for each clade using PhastSim on those topologies.
There are additional parameters near the top of the script (these will likely turn into command line parameters soon).

    num_res=x
    num_sim=y

These control randomness in how the simulated tree is generated.
`num_res` is the number of seeds to use for resolving multifurcations, and `num_sim` is the number of seeds to use for PhastSim.
The script will create `x * y` total "trials", each representing a simulation based on the given clade.
These simulations can be found at `data/<clade-name>/<trial-number>/simulation`.
In this directory, 
`ctree_with_refseq.fasta` contains the leaves of the simulated tree and the ancestral sequence and `collapsed_simulated_tree.nwk` contains the simulated tree.

### Inference

After running simulations, you can infer supports using 

    bash support_pipeline_scripts/run_inference_cluster.sh

This has similar parameters to the simulation script near the top of the file

    num_res=x
    num_sim=y

with the addition of

    method="beast"

which controls the method of inference. Currently, the script only supports `hdag` and `beast`.
This will produce a `results.pkl` in the directory

    data/<clade-name>/<trial-number>/results/<method>

for each clade and trial. This is a pickled list of tuples `(clade, estimate-support, in-tree)`, where

1. `clade` is the set of all taxa under a particular node.
2. `estimate-support` is the inferred support for that node.
3. `in-tree` is a boolean that indicates whether that node is in the true tree.

This direcotry also includes other information about the inference such as a log of trees for beast and the optimized dag.
