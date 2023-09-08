import sys
sys.path.append("/fh/fast/matsen_e/whowards/hdag-benchmark/support_pipeline")

import historydag as hdag
from historydag.utils import load_fasta, count_labeled_binary_topologies
from historydag.parsimony_utils import (
    compact_genome_hamming_distance_countfuncs
)
from historydag.compact_genome import (
    cg_diff
)
from historydag.utils import collapse_ete, resolve_ete, _test_resolve_ete
from historydag.parsimony import build_tree
from historydag.mutation_annotated_dag import CGHistoryDag
import pickle
import ete3
import time
from collections import Counter
from hdb import collapse_tree as ct
import numpy as np
import math
from scipy.linalg import expm
import os
import subprocess


experiment_type = "sim_models"
clade = sys.argv[1]
trial = sys.argv[2]
sim_type = "gamma_10_hmut_50" #sys.argv[3]
n_samples = 1000


dir_path = f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/{experiment_type}/{clade}/{sim_type}/{trial}"

fasta_path = dir_path + "/simulation/collapsed_simulated_tree.nwk.fasta"
dag_path = dir_path + "/results/historydag/final_opt_dag_trimmed.pb"
mut_rate_path = f"{dir_path}/simulation/sars-cov-2_simulation_output.info"

fasta = load_fasta(fasta_path)

def main():
    p_funcs = {
        "mut_rates_bl=1e-3": lambda n: prob_muts(n, branch_len=1e-3),
        "mut_rates_bl=1e-2": lambda n: prob_muts(n, branch_len=1e-2),
        "len_muts_scale=1e-3": lambda n: prob_muts_basic(n, scale=1e-3),
        "len_muts_scale=1e-2": lambda n: prob_muts_basic(n, scale=1e-2),
    }
    distr_dict = get_prob_muts_distribution(dag_path, p_funcs)
    with open(dir_path + f"/results/diff_historydag/prob_distrs.pkl", "wb") as f:
        pickle.dump(distr_dict, f)

                
def get_prob_muts_distribution(input_path, p_funcs):
    print("Loading and collpasing DAG...")
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(input_path)
    dag.convert_to_collapsed()
    dag.count_histories(bifurcating=True)
    dag.probability_annotate(
        edge_weight_func=lambda par, child: count_labeled_binary_topologies(len(child.clades))
    )
    seq2id = {seq: id for id, seq in fasta.items()}
    node2id = {}
    for node in dag.postorder():
        if node.is_leaf():
            seq = node.label.compact_genome.to_sequence()
            if seq in seq2id:
                node2id[node] = seq2id[seq]
            else:
                node2id[node] = "unnamed"
        else:
            node2id[node] = "unnamed"

    probs = {name: [] for name, _ in p_funcs.items()}
    for i, tree in enumerate(range(n_samples)):
        if i % (n_samples // 10) == 0:
            print(i)
        for name, p_func in p_funcs.items():
            history = dag.sample()
            tree = history.to_ete(
                name_func=lambda n: node2id[n], features=["compact_genome"]
            )
            for node in tree.traverse():
                if not node.is_root():
                    probs[name].append(p_func(node))
    return probs


_pb_nuc_lookup = {0: "A", 1: "C", 2: "G", 3: "T"}
_pb_nuc_codes = {nuc: code for code, nuc in _pb_nuc_lookup.items()}

substitution_model = np.array(
            [
                [-0.472, 0.039, 0.31, 0.123],
                [0.14, -3.19, 0.022, 3.028],
                [0.747, 0.113, -3.813, 2.953],
                [0.056, 0.261, 0.036, -0.353]
            ]
        )
site2rates = {}
with open(mut_rate_path, "r") as fh:
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

def prob_muts(ete_node, branch_len=0.0001):
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

def prob_muts_basic(ete_node, scale=0.01):
    muts = list(
                    cg_diff(
                        ete_node.up.compact_genome, ete_node.compact_genome
                    )
                )
    return scale * 1 / (len(muts)+1)


if __name__ == '__main__':
    main()