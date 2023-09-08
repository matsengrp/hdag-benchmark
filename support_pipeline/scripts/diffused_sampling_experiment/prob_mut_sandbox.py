import historydag as hdag
from historydag.parsimony import load_fasta, build_tree, disambiguate, sankoff_upward
from historydag.compact_genome import cg_diff
from historydag.utils import collapse_ete

import sys
import pickle
from collections import Counter

import numpy as np
from scipy.linalg import expm

import ete3

"""
Run:
bash support_pipeline/scripts/diffused_sampling_experiment/cluster_driver.sh prob_mut_sandbox.py
"""

experiment_type = "sim_models"
clade = sys.argv[1]
trial = sys.argv[2]
sim_type = "gamma_10_hmut_50"

_pb_nuc_lookup = {0: "A", 1: "C", 2: "G", 3: "T"}
_pb_nuc_codes = {nuc: code for code, nuc in _pb_nuc_lookup.items()}
dir_path = f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/{experiment_type}/{clade}/{sim_type}/{trial}"
mut_rates_path = f"{dir_path}/simulation/sars-cov-2_simulation_output.info"
sim_tree_path = f"{dir_path}/simulation/collapsed_simulated_tree.nwk"

fasta_path = dir_path + "/simulation/collapsed_simulated_tree.nwk.fasta"
dag_path = dir_path + "/results/historydag/final_opt_dag_trimmed.pb"
outdir = dir_path + "/results/diff_historydag_mut_rates"
fasta = load_fasta(fasta_path)

def get_rate_dict():
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
    return site2rates


def get_sim_tree_pcm_counts():
    """
    Returns two dictionarys mapping mutations to number of times it was seen in a PCM and the
    total number of times it was seen in the simulated tree.
    """
    tree = ete3.Tree(sim_tree_path, format=0)    # Load simulated tree with mutations on edges
    print(f"Tree contains {len(list(tree.get_descendants()))} nodes BEFORE collapse")
    tree = get_collapsed_sim_tree()
    print(f"Tree contains {len(list(tree.get_descendants()))} nodes AFTER collapse w/ hDAG")

    # Populate dictionary (from, to, site) -> (num_pcm, num_appeared)
    num_pcm = {}
    for node in tree.traverse():
        if node.is_leaf():
            continue
        
        # Count mutations on each child branch
        mut_counts = Counter()
        for child in node.children:
            muts = set(cg_diff(node.compact_genome, child.compact_genome))
            for mut in muts:
                mut_counts[mut] += 1

        # 
        for mut, count in mut_counts.items():
            if mut not in num_pcm:
                num_pcm[mut] = (0, 0)
            num_pcm[mut] = (num_pcm[mut][0], num_pcm[mut][1]+count)
            if count > 1:
                print("PCM FOUND:", mut, count)
                num_pcm[mut] = (num_pcm[mut][0]+count, num_pcm[mut][1])
    
    mut2num_pcm = {mut: count_pcm for mut, (count_pcm, count_observed) in num_pcm.items()}
    mut2count = {mut: count_observed for mut, (count_pcm, count_observed) in num_pcm.items()}
    return mut2num_pcm, mut2count

def get_collapsed_sim_tree():
    """ Loads the simulated tree for the given clade trial as a collapsed tree.
    """
    fasta = hdag.parsimony.load_fasta(f'{dir_path}/simulation/ctree_with_refseq.fasta')
    etetree = ete3.Tree(f'{dir_path}/simulation/collapsed_simulated_tree.nwk', format=0)

    def replace_site(site, newbase, sequence, oldbase=None):
        idx = site - 1
        if oldbase is not None and sequence[idx] != oldbase:
            print(f"WARNING: previous base was not as expected (expected={oldbase} vs actual={sequence[idx]} at {site})")
        return sequence[:idx] + newbase + sequence[idx + 1:]

    # Reconstruct sequences on each node
    ancestral_seq = fasta["ancestral"]
    for node in etetree.traverse("preorder"):
        if node.is_root():
            assert len(node.mutations) == 0
            parent_seq = ancestral_seq
        else:
            parent_seq = node.up.sequence
        if len(node.mutations) > 0:
            mut_list = node.mutations.split("|")
            for mut_str in mut_list:
                parent_seq = replace_site(int(mut_str[1:-1]), mut_str[-1], parent_seq, oldbase=mut_str[0])
        node.add_feature("sequence", parent_seq)

    seqs = set(fasta.values())
    for leaf in etetree.get_leaves():
        assert leaf.sequence in seqs

    sim_history = hdag.history_dag_from_trees(
        [etetree],
        [],
        label_functions = {"sequence": lambda n: n.sequence},
        attr_func = lambda n: {"name": n.name if n.is_leaf() else "internal"},
    )

    opt_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(dag_path)
    ref_seq = next(opt_dag.preorder(skip_ua_node=True)).label.compact_genome.reference
    sim_history = hdag.mutation_annotated_dag.CGHistoryDag.from_history_dag(sim_history, reference=ref_seq)
    sim_history.convert_to_collapsed()
    
    seq2id = {seq: id for id, seq in fasta.items()}
    node2id = {}
    for node in sim_history.postorder():
        if node.is_leaf():
            seq = node.label.compact_genome.to_sequence()
            if seq in seq2id:
                node2id[node] = seq2id[seq]
            else:
                node2id[node] = "unnamed"
        else:
            node2id[node] = "unnamed"

    ctree = sim_history.sample().to_ete(
                name_func=lambda n: node2id[n], features=["compact_genome"]
            )
    return ctree


def get_mut_counts():
    """
    Loads the optimized, MP-trimmed hdagReturns and returns a dictionary mapping mutations
    (from, to, site) -> 'num trees mutation appears in' / 'num trees'
    """
    opt_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(dag_path)
    mut_freq = {}  # (sequence_index, parent_nuc, child_nuc) -> frequency
    edge_counts = opt_dag.count_edges()    # Number of trees each edge takes part in
    total_muts = 0
    for child in reversed(list(opt_dag.postorder())):
        if not child.is_root():
            for parent in child.parents:
                if parent.is_root() or child.is_leaf():
                    continue
                muts = list(
                    cg_diff(
                        parent.label.compact_genome, child.label.compact_genome
                    )
                )

                for mut in muts:
                    if mut not in mut_freq:
                        mut_freq[mut] = 0
                    mut_freq[mut] += edge_counts[(parent, child)]
                    total_muts += edge_counts[(parent, child)]
    return mut_freq


def main():
    mut2num_pcm, mut2count = get_sim_tree_pcm_counts()
    print(f"A total of {len(mut2num_pcm)} types of mutations occured,",\
          f"{len([mut for mut, num_pcm in mut2num_pcm.items() if num_pcm > 0])} of which are present in a PCM")
    
    # site2rate = get_rate_dict()
    # prop_rate_list = []
    # for mut, num_pcm in mut2num_pcm.items():
    #     prop = num_pcm / mut2count[mut]
    #     from_nuc, to_nuc, site = mut
    #     rate_matrix = site2rate[site]
    #     from_idx = _pb_nuc_codes[from_nuc]
    #     to_idx = _pb_nuc_codes[to_nuc]
    #     mut_rate = rate_matrix[from_idx, to_idx]
    #     prop_rate_list.append((prop, mut_rate))
    # with open(outdir + "/num_pcm_rate_list.pkl", "wb") as f:
    #     pickle.dump(prop_rate_list, f)

    # mut_freqs = get_mut_counts()
    # print(list(mut_freqs.keys()))
    # prop_freq_list = []
    # for mut, prop in prop_pcm.items():
    #     if mut not in mut_freqs:
    #         mut_freq = 0
    #     else:
    #         mut_freq = mut_freqs[mut]
    #     prop_freq_list.append((prop, mut_freq))
    # with open(outdir + "/num_pcm_freq_list.pkl", "wb") as f:
    #     pickle.dump(prop_freq_list, f)


    site2rate = get_rate_dict()
    counts_rate_list = []
    for mut, num_pcm in mut2num_pcm.items():
        count = mut2count[mut]

        # Get the site for the given mutation
        from_nuc, to_nuc, site = mut
        rate_matrix = site2rate[site]
        from_idx = _pb_nuc_codes[from_nuc]
        to_idx = _pb_nuc_codes[to_nuc]
        mut_rate = rate_matrix[from_idx, to_idx]

        counts_rate_list.append((num_pcm, count, mut_rate))
    
    with open(outdir + "/counts_rate_list.pkl", "wb") as f:
        pickle.dump(counts_rate_list, f)

    
    






if __name__ == '__main__':
    main()