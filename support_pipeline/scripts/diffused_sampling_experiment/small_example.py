import os
import sys
import subprocess
import pickle
sys.path.append("/fh/fast/matsen_e/whowards/hdag-benchmark/support_pipeline")
from utils import reroot

import ete3
import pydot
import numpy as np
from sklearn.manifold import TSNE, MDS

import historydag as hdag
from historydag.mutation_annotated_dag import CGHistoryDag
from collections import Counter

from historydag.utils import (
    collapse_ete,
    resolve_ete,
    count_labeled_binary_topologies,
)

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
mpl.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 17,
    'axes.labelsize': 16,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'font.family': 'sans-serif',
    'font.weight': 600,
    'axes.labelweight': 600,
    'axes.titleweight': 600,
    'figure.autolayout': True
    })

import random
random.seed(123)


local_outpath = "/fh/fast/matsen_e/whowards/hdag-benchmark/tmp"

def small_example():
    """
    Generates a small example dataset.
    """

    fasta = {
        "ancestral": "AAAAA",
        "s1": "GTCCC",
        "s2": "GTCAT",
        "s3": "GACAG",
        "s4": "ATCAA",
        "s5": "ATAAC",
    }
    seq2id = {seq: id for id, seq in fasta.items()}
    tree = ete3.Tree()
    tree.add_feature("sequence", fasta["ancestral"])
    tree.name = "ancestral"
    for name, seq in fasta.items():
        if name != "ancestral":
            child = ete3.Tree()
            child.add_feature("sequence", seq)
            child.name = name
            tree.add_child(child)

    opt_dag = get_opt_dag(tree, fasta)
    dist_mat, bifurcating_topologies, history_list = get_dist_mat(tree, fasta, opt_dag)
    bif_node2id = get_node2id(bifurcating_topologies, fasta)
    node2id = get_node2id(opt_dag, fasta)
    ete_tree_list = [history.to_ete(name_func=lambda n: bif_node2id[n]) for history in history_list]

    n_samples = 10000

    counts_list = []
    props = [0, 0.05, 0.1, 0.15, 0.2, 0.25]
    # for prop in props:
    #     diff_tree_counts = Counter()
    #     diffused_sampler = opt_dag.copy().diffused_tree_sampler(n_samples, node2id, lambda n: prop)
    #     for i in range(n_samples):
    #         if i % int(n_samples/10) == 0:
    #             print("sample:", prop, i)
    #         rtree = next(diffused_sampler)

    #         dist = np.zeros(len(ete_tree_list))
    #         for j, tree in enumerate(ete_tree_list):
    #             dist[j] = rf_dist(rtree, tree)
    #         idx = np.where((dist_mat == dist).all(axis=1))[0][0]
    #         diff_tree_counts[idx] += 1
    #     counts_list.append(diff_tree_counts)

    #     print(f"=> Sampled {len(diff_tree_counts)} unique trees from the {prop} distribution")
    
    # with open(f"{local_outpath}/counts_list.pkl", "wb") as fh:
    #     pickle.dump(counts_list, fh)

    with open(f"{local_outpath}/counts_list.pkl", "rb") as fh:
        counts_list = pickle.load(fh)


    print("Plotting embeddings...")
    plot_embs(dist_mat, counts_list, name="all", titles=props)

    
    # TODO: Figure out why the MP distribution is so centered around the tree at index 58
    print("=> Examining the \"uniform\" distirbution on MP trees")
    prop = props[0]
    counts = counts_list[0]
    print(len(counts), "unique trees sampled. Here are their counts:")
    max_count = -1
    max_idx = -1
    for idx, count in counts.items():
        if count > max_count:
            max_count = count
            max_idx = idx
        print(idx, count)
    print()

    # NOTE: Checking resolve_ete
    # total_trees = 0
    # trials = 10000
    # opt_dag.convert_to_collapsed()
    # for i, history in enumerate(opt_dag.get_histories()):
    #     ctree = history.to_ete(name_func=lambda n: node2id[n], features=["compact_genome"])
    #     print(f"{i+1}. RF dist to tree #{max_idx}", rf_dist(ete_tree_list[max_idx], ctree))
    #     tree_list = []
    #     tree_counts = Counter()
    #     num_toi_samples = 0
    #     for t in range(trials):
    #         rtree = ctree.copy()
    #         resolve_ete(rtree)
    #         idx = len(tree_list)
    #         for j, tree in enumerate(tree_list):
    #             if ete_equals(tree, rtree):
    #                 idx = j
    #                 break
    #         if idx == len(tree_list):
    #             tree_list.append(rtree)
    #         tree_counts[idx] += 1
    #         if rf_dist(ete_tree_list[max_idx], rtree) == 0:
    #             num_toi_samples += 1
    #     print(f"\t{len(tree_counts)} unique trees. ToI sampled {num_toi_samples}/{trials} times")
    #     total_trees += len(tree_counts)
    #     for tree in tree_list:
    #         print(f"\tRF dist to tree #{max_idx}", rf_dist(ete_tree_list[max_idx], tree))
    #         # print(tree)
    #     print()
    #     print("Total Trees:", total_trees)

    # TODO: Check the proportion from hdag.sample()
    trials = 1000
    opt_dag.convert_to_collapsed()
    num_trees = opt_dag.count_histories(bifurcating=True)

    # TODO: Why does this put so much probability mass on tree #0?
    opt_dag.probability_annotate(
        edge_weight_func=lambda par, child: count_labeled_binary_topologies(len(child.clades))
    )

    tree_list = []
    tree_counts = Counter()
    for i in range(trials):
        history = opt_dag.sample()
        ctree = history.to_ete(name_func=lambda n: node2id[n], features=["compact_genome"])
        # print(f"{i+1}. RF dist to tree #{max_idx}", rf_dist(ete_tree_list[max_idx], ctree))
        idx = len(tree_list)
        for j, tree in enumerate(tree_list):
            if ete_equals(tree, ctree):
                idx = j
                break
        if idx == len(tree_list):
            tree_list.append(ctree)
        tree_counts[idx] += 1

    print(f"DAG contains {opt_dag.count_histories()} histories ({num_trees} bifurcating)")
    print(f"\t{len(tree_counts)} unique trees.")
    for i, tree in enumerate(tree_list):
        num_resolutions = 1
        for node in tree.traverse():
            num_resolutions *= count_labeled_binary_topologies(len(node.children))
        print(f"{i}. Has {num_resolutions} resolutions, but observed in {tree_counts[i]}/{trials} samples", tree)
    print()
    print(tree_counts)


def plot_embs(dist_matrix, list_idx2count, name="", titles=None):
    fig, axes = plt.subplots(nrows=1, ncols=5, sharey=True, sharex=True)

    embedding = MDS(n_components=2, dissimilarity='precomputed', random_state=123)
    # embedding = TSNE(n_components=2, metric="precomputed", init='random', random_state=123)
    coords = embedding.fit_transform(dist_matrix)
    x = np.array([coord[0] for coord in coords])
    y = np.array([coord[1] for coord in coords])

    cmap = cm.get_cmap("viridis")
    fig.set_size_inches(w=18, h=4)

    max_prob = max(list_idx2count[0].values())
    total = sum(list_idx2count[0].values())
    max_prob /= total

    print("\tMax prob:", max_prob)

    dag_idxs = None
    for i, (idx2count, ax) in enumerate(zip(list_idx2count, axes)):
        prob = np.zeros(len(dist_matrix))
        for j, count in idx2count.items():
            prob[j] = count
        prob /= total
        sampled_idxs = np.nonzero(prob)[0]
        if dag_idxs is None:
            dag_idxs = sampled_idxs
        sampled_idxs = np.array([i for i in sampled_idxs if i not in list(dag_idxs)])
        if len(sampled_idxs) > 0:
            ax.scatter(x[sampled_idxs], y[sampled_idxs], c=prob[sampled_idxs], cmap=cmap, vmin=0, vmax=max_prob, s=50, alpha=0.75)
        ax.scatter(x[dag_idxs], y[dag_idxs], c=prob[dag_idxs], cmap=cmap, vmin=0, vmax=max_prob, marker='*', s=100)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        if titles is not None:
            ax.set_title(titles[i])

    # inches = 1.5
    # plt.subplots_adjust(right=inches)
    # cax = fig.add_axes([inches + 0.05, 0.15, 0.05, 0.7])
    # plt.colorbar(im, cax=cax)

    # fig.colorbar(im, ax=axes.ravel().tolist())

    # cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
    # plt.colorbar(im, cax=cax, **kw)

    fig.suptitle("t-SNE Projection of Pairwise RF Distances")
    plt.tick_params(left = False, bottom = False)
    plt.savefig(f"{local_outpath}/tree_embs_{name}.pdf")
    plt.clf()


def get_opt_dag(tree, fasta):
    """
    Builds a DAG with MP trees on the given fasta.
    """
    # with open(local_outpath + "/var_sites_with_refseq.fasta", "w") as f:
    #     for name, seq in fasta.items():
    #         txt = f">{name}\n{seq}\n"
    #         f.write(txt)
    # history = CGHistoryDag.from_history_dag(hdag.history_dag_from_etes([tree], ['sequence']), reference=fasta["ancestral"])
    # larch(history, local_outpath, 2000, pscore_coeff=1, node_coeff=2)
    # opt_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(str(local_outpath) + "/trimmed_opt_dag.pb")
    # print(opt_dag.count_histories())

    # subprocess.run([
    #     "bash",
    #     "/fh/fast/matsen_e/whowards/hdag-benchmark/hdb/scripts/dnapars_parsimony_score.sh",
    #     local_outpath + "/var_sites_with_refseq.fasta",          # infasta
    #     "ancestral",                                                   # root_name
    #     local_outpath                                           # out_dir
    # ])

    # with open(local_outpath + "/dnapars_output.txt", "r") as f:
    #     line = f.readline()
    #     temp = line.strip().split(" ")
    #     best_possible = float(temp[-1])
    #     print("Best pscore:", best_possible)

    # dnapars_tree_list = []
    # with open(local_outpath + "/dnapars/outtree", "r") as f:
    #     for line in f:
    #         if line[-1] != ';':
    #             nwk = line.strip() + next(f).strip().split("[")[0] + ";"
    #         tree = ete3.Tree(nwk, format=5)
    #         rerooted = reroot(tree.search_nodes(name="ancestral")[0])
    #         for node in rerooted.traverse():
    #             if node.is_root():
    #                 node.name = "ancestral"
    #                 node.add_feature("sequence", fasta["ancestral"])
    #             elif node.is_leaf():
    #                 node.add_feature("sequence", fasta["s" + node.name[3:]])
    #             else:
    #                 # node.name = "**"
    #                 node.add_feature("sequence", fasta["ancestral"])
    #         dnapars_tree_list.append(rerooted)

    # temp_dag = hdag.history_dag_from_etes(dnapars_tree_list, ['sequence'])
    # dnapars_dag = CGHistoryDag.from_history_dag(temp_dag, reference=fasta["ancestral"])
    # opt_dag.merge(dnapars_dag)

    # print("Before Larch opt there are:", opt_dag.count_histories())

    # node2id = {}
    # seq2id = {v:k for k, v in fasta.items()}
    # for node in opt_dag.postorder():
    #     if node.is_leaf():
    #         seq = node.label.compact_genome.to_sequence()
    #         node2id[node] = seq2id[seq]
    #     else:
    #         node2id[node] = "unnamed"
    # for history in opt_dag.get_histories():
    #     tree = history.to_ete(name_func=lambda n: node2id[n], features=["compact_genome"])
    #     print(tree)

    # larch(opt_dag, local_outpath, 100, pscore_coeff=1, node_coeff=2)  # TODO: This errors out after >100 iters
    opt_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(str(local_outpath) + "/trimmed_opt_dag.pb")
    opt_dag.make_complete()
    opt_dag.convert_to_collapsed()
    opt_dag.trim_optimal_weight()
    print("After Larch opt there are:", opt_dag.count_histories())
    return opt_dag

def get_dist_mat(tree, fasta, dag=None):
    """
    Generate the distance matrix on all bifurcating topologies on the taxa that `tree` is on.
    """
    tree_list = []
    n_leaves = len(list(tree.get_leaves()))
    n_samples = 10 * count_labeled_binary_topologies(n_leaves)
    for i in range(n_samples):
        rtree = tree.copy()
        resolve_ete(rtree, feature_func_dict={"sequence": lambda n: n.sequence})
        tree_list.append(rtree)
    bifurcating_histories = CGHistoryDag.from_history_dag(hdag.history_dag_from_etes(tree_list, ['sequence']), reference=fasta["ancestral"])
    
    if dag is not None:
        cg_tree_list = []
        node2id = get_node2id(dag, fasta)
        mp_sampler = dag.diffused_tree_sampler(n_samples, node2id, lambda n: 0)
        for i in range(1000):
            rtree = next(mp_sampler)
            for node in rtree.traverse():
                assert node.is_leaf() or len(node.children) == 2
            # cg_history = CGHistoryDag.from_history_dag(history, reference=fasta["ancestral"])
            cg_tree_list.append(rtree)
        cg_histories = hdag.history_dag_from_etes(cg_tree_list, ['compact_genome'])
        bifurcating_histories.merge(cg_histories)

    bifurcating_histories = bifurcating_histories.unlabel()
    bifurcating_histories.make_complete()

    num_top = bifurcating_histories.count_histories()
    num_leaves = 0
    for node in bifurcating_histories.postorder():
        if node.is_leaf():
            num_leaves += 1

    print(f"After sampling {n_samples} times there are {num_top} topologies in our DAG on {num_leaves} leaves")
    
    node2id = get_node2id(bifurcating_histories, fasta)
    dist = np.zeros((num_top, num_top))
    tree_list = [t for t in bifurcating_histories]
    for i, t1 in enumerate(tree_list):
        if i % 10 == 0:
            print(i)
        rf_dist_func = hdag.utils.make_rfdistance_countfuncs(t1)
        for j, t2 in enumerate(tree_list[i+1:], i+1):

            # TODO: Figure out which RF distance is correct
            # rf = t2.optimal_weight_annotate(**rf_dist_func)
            
            t1_ete = t1.to_ete(name_func=lambda n: node2id[n])
            t2_ete = t2.to_ete(name_func=lambda n: node2id[n])
            ete_rf = rf_dist(t1_ete, t2_ete)
            
            # if ete_rf != rf:
            #     print(i, j, "\t", ete_rf, "\t", rf)
            #     # IDEA: Print these trees and compare manaully
            # if i == 59 and j == 78:
            #     print("t1 vs t2")
            #     print(t1_ete)
            #     print(t2_ete)

            dist[i, j] = ete_rf
            dist[j, i] = ete_rf

    return dist, bifurcating_histories, tree_list

def get_node2id(dag, fasta):
    # NOTE: Only works on CG labelled DAG
    node2id = {}
    seq2id = {v:k for k, v in fasta.items()}
    for node in dag.postorder():
        if node.is_leaf():
            seq = node.label.compact_genome.to_sequence()
            node2id[node] = seq2id[seq]
        else:
            node2id[node] = "unnamed"
    return node2id

    
def ete_equals(t1, t2):
    """Uses RF distance to determine whether two ete trees are equal"""
    return rf_dist(t1, t2) == 0

def rf_dist(t1, t2):
    assert set([l.name for l in t1.get_leaves()]) == set([l.name for l in t2.get_leaves()])
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

    rf_dist = len(set(t1_node2cu.values()).difference(set(t2_node2cu.values()))) + \
            len(set(t2_node2cu.values()).difference(set(t1_node2cu.values())))
    return rf_dist

def larch(input, out_dir, count=1000, node_coeff=3, pscore_coeff=1, executable=None):
    """Python method for driving larch"""
    if executable is None:
        executable = '/home/whowards/larch/larch/build/larch-usher'

    if isinstance(input, str):
        print("Using string as input")
    else:
        input_path = f"{out_dir}/sim_dag.pb"
        print(f"Assuming input is historyDAG. Saving to disk at {input_path}...")
        input.to_protobuf_file(input_path)
        input = input_path

    os.chdir(f"{out_dir}") # E.g., results/historydag/
    log_dir = f"{out_dir}/opt_info/optimization_log"

    print(f"Running {count} iterations of Larch...")
    subprocess.run(["mkdir", "-p", f"{log_dir}"])
    args = [executable,
            "-i", f"{input}",
            "-c", f"{count}",
            "-o", f"{out_dir}/trimmed_opt_dag.pb",
            "-l", f"{log_dir}",
            "--move-coeff-nodes", str(node_coeff),
            "--move-coeff-pscore", str(pscore_coeff),
            "--sample-any-tree",
            ]
    with open(f"{log_dir}/mat_opt.log", "w") as f:
        subprocess.run(args=args, stdout=f, stderr=f)



def main():
    small_example()

if __name__ == '__main__':
    main()



