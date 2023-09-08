import sys
sys.path.append("/fh/fast/matsen_e/whowards/hdag-benchmark/support_pipeline")

from utils import reroot, larch

import historydag as hdag
from historydag.utils import load_fasta
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

n_samples = 100

base_path = "/fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models/AY.74/gamma_10_hmut_50/7"
inp = f"{base_path}/results/historydag/final_opt_dag_trimmed.pb"
id2seq = load_fasta(f"{base_path}/simulation/ctree_with_refseq.fasta")
# input_newick = f"{base_path}/simulation/sars-cov-2_simulation_output.tree"
mut_rate_path = f"{base_path}/simulation/sars-cov-2_simulation_output.info"
ctree_path = f"{base_path}/simulation/collapsed_simulated_tree.nwk"
mb_tree_file = f"{base_path}/results/mrbayes/mrbayes-output.t"

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
site2lnrate = {}
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
        
        site2lnrate[site] = Q

def prob_muts(ete_node):
    muts = list(
                    cg_diff(
                        ete_node.up.compact_genome, ete_node.compact_genome
                    )
                )
    prob_muts = 1
    t = 0.001 
    for old, new, site in muts:
        i = _pb_nuc_codes[old]
        j = _pb_nuc_codes[new]
        Q = site2lnrate[site]
        assert all([np.isclose(np.sum(Q[i]), 0.0) for i in range(4)])
        prob_mut = expm(t*Q)[i, j]
        assert prob_mut <= 1 and prob_mut >= 0
        prob_muts *= prob_mut
    return prob_muts

def prob_muts_basic(ete_node):
    muts = list(
                    cg_diff(
                        ete_node.up.compact_genome, ete_node.compact_genome
                    )
                )
    return 0.01 * 1 / len(muts)

def rf_dist(t1, t2):
    t1_node2cu = {}
    for node in t1.traverse("postorder"):
        if node.is_leaf():
            cu = [node.name]
        else:
            cu = []
            for child in node.children:
                cu.extend(list(t1_node2cu[child]))
        cu = frozenset(cu)
        t1_node2cu[node] = cu
    t2_node2cu = {}
    for node in t2.traverse("postorder"):
        if node.is_leaf():
            cu = [node.name]
        else:
            cu = []
            for child in node.children:
                cu.extend(list(t2_node2cu[child]))
        cu = frozenset(cu)
        t2_node2cu[node] = cu
    return len(set(t1_node2cu.values()).difference(set(t2_node2cu.values()))) + \
            len(set(t2_node2cu.values()).difference(set(t1_node2cu.values())))


def get_mrbayes_trees(trees_file, skip=0):
    with open(trees_file, 'r') as fh:
        for line in fh:
            if 'translate' in line:
                break
        translate_dict = {}
        for line in fh:
            idx, idx_name = line.strip().split(' ')
            translate_dict[idx] = idx_name[:-1]
            if idx_name[-1] == ';':
                break

    with open(trees_file, 'r') as fh:
        for i, line in enumerate(fh):
            if skip > 0 and i % skip:
                continue
            if 'tree' in line.strip()[0:5]:
                nwk = line.strip().split(' ')[-1]
                tree = build_tree(nwk, translate_dict)
                for node in tree.iter_leaves():
                    node.name = translate_dict[node.name]
                rerooted = reroot(tree.search_nodes(name="ancestral")[0])
                yield rerooted

# def get_trees(input_path):
#     i = -1
#     leaf_map = {}
#     with open(input_path, 'r') as fh:
#         for line in fh:
#             if '    ' in line[0:5]:
#                 mb_leaf, fasta_leaf = tuple(line.strip().split(' '))
#                 leaf_map[mb_leaf] = fasta_leaf[:-1]
#                 fl2ml = {v: k for k, v in leaf_map.items()}
#             elif 'tree' in line.strip()[0:5]:
#                 i += 1
#                 if i < burnin * num_trees:
#                     continue
#                 treeprob = 1
#                 nwk = line.strip().split(' ')[-1]

#                 fasta_mb = {fl2ml[fasta_leaf]: v for fasta_leaf, v in fasta.items()}
#                 tree = build_tree(nwk, fasta_mb)

#                 pars_score = sankoff_upward(
#                     tree,
#                     len(tree.sequence),
#                     transition_model=hdag.parsimony_utils.default_nt_transitions,
#                     use_internal_node_sequences=False,
#                 )
#                 tree.pars_score = pars_score
#                 tree.prob = treeprob

#                 if i % 10 == 0:
#                     print(i, pars_score)

#                 yield tree

def main():
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(inp)
    dag.convert_to_collapsed()
    with open(ctree_path, 'r') as fh:
        true_tree = ete3.Tree(fh.read())

    # TODO: Compute best possible parsimony score for this topology?
    n_leaves = len(true_tree.get_leaves())
    true_pscore = 0
    true_num_non_mut_nodes = 0
    for node in true_tree.traverse():
        true_pscore += len(node.mutations)
        true_num_non_mut_nodes += int(len(node.mutations) == 0)

    print(f"True tree has a parsimony score of {true_pscore}")


    seq2id = {seq: id for id, seq in id2seq.items()}
    node2id = {}
    for node in dag.postorder():
        if node.is_leaf():
            seq = node.label.compact_genome.to_sequence()
            if seq in seq2id:
                node2id[node] = seq2id[seq]
            else:
                node2id[node] = "unnamed"

    scores = dag.hamming_parsimony_count()
    pscore = next(iter(scores.keys()))
    print("Sampling from MP-Diffused distribution with MP =", pscore)
    print(f"True tree has pscore of {true_pscore} and {true_num_non_mut_nodes} / {2 * n_leaves - 2} branches w/o mutations")

    avg_mb_dist = 0
    avg_opt_mp_dist = 0
    avg_opt_diff_dist = 0
    avg_mp_dist = 0
    avg_diff_dist = 0
    avg_rand_dist = 0
    mp_dag = dag.copy()

    sampler = dag.diffused_tree_sampler(n_samples, node2id, prob_muts)

    optimal_dag = get_sim_tree_optimized_dag()

    opt_mp_dag = optimal_dag.copy()
    opt_node2id = {}
    for node in opt_mp_dag.postorder():
        if node.is_leaf():
            seq = node.label.compact_genome.to_sequence()
            if seq in seq2id:
                opt_node2id[node] = seq2id[seq]
            else:
                opt_node2id[node] = "unnamed"
        else:
            opt_node2id[node] = "unnamed"

    opt_sampler = optimal_dag.diffused_tree_sampler(n_samples, opt_node2id, prob_muts_basic)

    num_mb_trees = 100000
    skip = 100
    burnin = int(0.5 * num_mb_trees / skip)
    mb_sampler = get_mrbayes_trees(mb_tree_file, skip=skip)
    for i in range(burnin):
        next(mb_sampler)

    for i in range(n_samples):
        opt_mp_history = opt_mp_dag.sample()
        opt_mp_tree = opt_mp_history.to_ete(name_func=lambda n: opt_node2id[n], features=["compact_genome"])
        resolve_ete(opt_mp_tree)

        opt_tree = next(opt_sampler)
        opt_sample_history = hdag.history_dag_from_etes([opt_tree], ['compact_genome'])
        opt_scores = opt_sample_history.weight_count(**compact_genome_hamming_distance_countfuncs)

        mp_history = mp_dag.sample()
        mp_tree = mp_history.to_ete(name_func=lambda n: node2id[n], features=["compact_genome"])
        resolve_ete(mp_tree)

        tree = next(sampler)
        sample_history = hdag.history_dag_from_etes([tree], ['compact_genome'])
        scores = sample_history.weight_count(**compact_genome_hamming_distance_countfuncs)
        
        random_tree = true_tree.copy()
        collapse_ete(random_tree, prob_muts=lambda n: 1)
        assert len(random_tree.children) == n_leaves
        resolve_ete(random_tree)
        assert len(random_tree.children) == 2

        mb_tree = next(mb_sampler)

        rf_mb_to_true = rf_dist(true_tree, mb_tree)
        rf_opt_mp_sampled_to_true = rf_dist(true_tree, opt_mp_tree)
        rf_opt_sampled_to_true = rf_dist(true_tree, opt_tree)
        rf_mp_sampled_to_true = rf_dist(true_tree, mp_tree)
        rf_sampled_to_true = rf_dist(true_tree, tree)
        rf_random_to_true = rf_dist(true_tree, random_tree)

        avg_mb_dist += rf_mb_to_true
        avg_opt_mp_dist += rf_opt_mp_sampled_to_true
        avg_opt_diff_dist += rf_opt_sampled_to_true
        avg_mp_dist += rf_mp_sampled_to_true
        avg_diff_dist += rf_sampled_to_true
        avg_rand_dist += rf_random_to_true
        print(f"{i}\tmb-dist: {rf_mb_to_true}\topt-diff-sample pscore: {next(iter(opt_scores.keys()))}\tdiff-sample pscore: {next(iter(scores.keys()))}\t opt-mp-sample rf:{rf_opt_mp_sampled_to_true}\t opt-diff-sample rf:{rf_opt_sampled_to_true}\t mp-sample rf:{rf_mp_sampled_to_true}\t diff-sample rf:{rf_sampled_to_true}\t random rf:{rf_random_to_true}")

    # print(f"Estimated probability of true tree...\n \
    #       => from OPT MP distribution is: {total_mp_sample/n_samples}\n \
    #       => from OPT MP-diffused distribution is: {total_sample/n_samples}\n \
    #       => from MP distribution is: {total_mp_sample/n_samples}\n \
    #       => from MP-diffused distribution is: {total_sample/n_samples}\n \
    #       => from unifrorm distribution is: {total_random/n_samples}")

    print(f"Average RF distance to true tree\n \
          => MrBayes: {avg_mb_dist / n_samples}\n \
          => OPT MP: {avg_opt_mp_dist / n_samples}\n \
          => OPT Diff: {avg_opt_diff_dist / n_samples}\n \
          => MP: {avg_mp_dist / n_samples}\n \
          => Diff: {avg_diff_dist / n_samples}\n \
          => Rand: {avg_rand_dist / n_samples}")
    
    print()



if __name__ == '__main__':
    main()
