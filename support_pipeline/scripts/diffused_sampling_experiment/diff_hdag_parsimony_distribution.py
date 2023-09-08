import historydag as hdag
from historydag.parsimony import load_fasta, build_tree, disambiguate, sankoff_upward
from historydag.compact_genome import cg_diff

import sys
import pickle
from collections import Counter

import numpy as np
from scipy.linalg import expm

"""
Run:
bash support_pipeline/scripts/diffused_sampling_experiment/cluster_driver.sh diff_hdag_parsimony_distribution
"""

experiment_type = "sim_models"
clade = sys.argv[1]
trial = sys.argv[2]
sim_type = "gamma_10_hmut_50" #sys.argv[3]

_pb_nuc_lookup = {0: "A", 1: "C", 2: "G", 3: "T"}
_pb_nuc_codes = {nuc: code for code, nuc in _pb_nuc_lookup.items()}
dir_path = f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/{experiment_type}/{clade}/{sim_type}/{trial}"
mut_rates_path = f"{dir_path}/simulation/sars-cov-2_simulation_output.info"

fasta_path = dir_path + "/simulation/collapsed_simulated_tree.nwk.fasta"
dag_path = dir_path + "/results/historydag/final_opt_dag_trimmed.pb"
outdir = dir_path + "/results/diff_historydag_mut_rates"
fasta = load_fasta(fasta_path)

def main():
    distr = get_parsimony_distribution(dag_path)
    print(distr)
    with open(outdir + f"/pars_distr.pkl", "wb") as f:
        pickle.dump(distr, f)

def prob_muts_basic(ete_node):
    muts = list(
                    cg_diff(
                        ete_node.up.compact_genome, ete_node.compact_genome
                    )
                )
    return 0.01 * 1 / (len(muts)+1)

def prob_muts(ete_node, site2rates, branch_len=1e-3):
    muts = list(
                    cg_diff(
                        ete_node.up.compact_genome, ete_node.compact_genome
                    )
                )
    prob_muts = 1

    for old, new, site in muts:
        i = _pb_nuc_codes[old]
        j = _pb_nuc_codes[new]

        Q = site2rates[site]
        
        assert all([np.isclose(np.sum(Q[i]), 0.0) for i in range(4)])
        prob_mut = expm(branch_len*Q)[i, j]
        assert prob_mut <= 1 and prob_mut >= 0
        prob_muts *= prob_mut
    
    return prob_muts
                
def get_parsimony_distribution(input_path):
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(input_path)
    dag.convert_to_collapsed()

    substitution_model = np.array(
            [
                [-0.472, 0.039, 0.31, 0.123],
                [0.14, -3.19, 0.022, 3.028],
                [0.747, 0.113, -3.813, 2.953],
                [0.056, 0.261, 0.036, -0.353]
            ]
        )
    site2rates = {}
    if mut_rates_path is not None:
        with open(mut_rates_path, "r") as fh:
            lines = fh.readlines()
            for line in lines[1:]:
                arr = line.strip().split("\t")

                # Adjust based on gamma rate variation
                site = int(arr[0])
                rate = float(arr[1])
                rate_mat = np.ones((4, 4)) * rate

                # Adjust based on hypermutation
                if arr[2] == "1":
                    i = _pb_nuc_codes[arr[3]]
                    j = _pb_nuc_codes[arr[4]]
                    rate_mat[i, j] *= 50

                Q = np.multiply(substitution_model, rate_mat)
                for i in range(4):
                    Q[i, i] -= np.sum(Q[i])
                
                site2rates[site] = Q

    n_samples = 10000
    branch_len = 0.001
    sampler = dag.diffused_tree_sampler(n_samples, fasta, lambda n: prob_muts(n, site2rates, branch_len=branch_len))

    distribution = Counter()
    for tree in sampler:
        pars = 0
        for node in tree.traverse():
            if not node.is_root():
                pars += len(list(cg_diff(
                        node.up.compact_genome, node.compact_genome
                    )))

        # TODO: Compute best parsimony score this:
        # pars_score = sankoff_upward(
        #             tree,
        #             len(tree.sequence),
        #             transition_model=hdag.parsimony_utils.default_nt_transitions,
        #             use_internal_node_sequences=False,
        #         )
        distribution[pars] += 1
    return distribution


if __name__ == '__main__':
    main()