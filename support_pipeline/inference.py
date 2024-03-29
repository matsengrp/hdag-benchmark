import click
import re
import ete3 as ete
import random
import pickle
import os

import random
import subprocess

# from newick_parser.tree_transformer import iter_nexus_trees
from historydag import parsimony_utils
from math import exp
import numpy as np
from scipy.linalg import expm

import time

import historydag as hdag
from historydag import parsimony_utils
from historydag.parsimony import build_tree
from historydag.utils import count_labeled_binary_topologies, load_fasta
from historydag.compact_genome import cg_diff

from utils import get_true_nodes, make_results_list, reroot, get_history

_pb_nuc_lookup = {0: "A", 1: "C", 2: "G", 3: "T"}
_pb_nuc_codes = {nuc: code for code, nuc in _pb_nuc_lookup.items()}

@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Command line scripts to facilitate node support validation
    """
    pass


@click.command("save_supports")
@click.option('--method', '-m', default="hdag", help='method of inference.')
@click.option('--tree_path', '-t', help='the newick file for the true tree.')
@click.option('--fasta_path', '-f', default="", help='the fasta file for the true tree (optional).')
@click.option('--input_path', '-i', help='the file containing inference results.')
@click.option('--output_path', '-o', help='the file to save output to.')
@click.option('--use_results', '-u', is_flag=True, help='edit the results file already stored at output path.')
def save_supports(method, tree_path, fasta_path, input_path, output_path, use_results):
    """
    A method for computing, formatting, and saving node supports for various methods of inference 
    Input:  method of inference (e.g., hdag, beast, dag-inf, hdag-adj, etc.),
            input file path (e.g., to history dag, or beast output),
            output filepath 
    Output: For each file in the input directory, a text file containing list of support values in the
            format (clade, estimated_support, in_tree)

    E.g.
        python support_pipeline_scripts/cli.py save_supports \
        -m hdag-inf \
        -t /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.132/2/simulation/collapsed_simulated_tree.nwk \
        -i /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.132/2/results/historydag/final_opt_dag.pb \
        -o /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.132/2/results/historydag/results_adj.pkl
    """

    # Map of taxon id (e.g., s1, s4, etc) to full sequence

    if fasta_path == "":
        fasta_path = tree_path + ".fasta"
    taxId2seq = hdag.parsimony.load_fasta(fasta_path)

    # Compute the set of nodes that are in the true tree
    node_set = get_true_nodes(tree_path, taxId2seq)

    if use_results:
        print("=> Using previous support results")
        with open(output_path, "rb") as f:
            results = pickle.load(f)
        support = {result[0]: result[1] for result in results}
        support_list = make_results_list(support, node_set)
    else:
        print("=> Computing new list of node supports")
        # Computes list of node supports
        if method == "hdag":
            p = "inf"
            adjust = False
            support_list = hdag_output_general(node_set, input_path, taxId2seq, pars_weight=p, adjust=adjust)#, true_tree_path=tree_path)
        elif method == "hdag-diff":
            support_list = nwk_sample_output(node_set, input_path)
        elif method[0:5] == "hdag-":
            if method[5:] == "inf":
                p = "inf"
                adjust = False
            elif method[5:] == "adj":
                p = "inf"
                adjust = True
            else:
                p = float(method[5:])
                adjust = False
            support_list = hdag_output_general(node_set, input_path, taxId2seq, pars_weight=p, adjust=adjust)
        elif method == "mrbayes":
            support_list = mrbayes_output(node_set, input_path)
        elif method == "ufboot":
            support_list = ufboot_output(node_set, input_path)
        elif method == "beast":
            raise Exception("BEAST inference is currently not supported.")
        # elif method == "random":
        #     support_list = random_output(node_set, input_path)
        else:
            raise Exception(f"Invalid method: {method}")

    # NOTE: Some
    # support_list = [result for result in support_list if result[1] > 0]
    print(len([result for result in support_list if result[1] == 0]), "nodes with 0 support")

    with open(output_path, "wb") as f:
        pickle.dump(support_list, f)
    
    # TODO: Debug the creation of the node set from the tree!
    #       Bug found! Was not changing name of resolved node
    # # But they aren't in the nodeset??
    # for node in node_set:
    #     print(node)

    # # Many modes with in_tree = False are actually in the true tree
    # no_leaves = [el for el in support_list if len(el[0]) > 1]
    # for el in no_leaves[-50:]:
    #     print(f"{el[2]}\t{el[0]}\t{el[1]}")


###################################################################################################
###     hDAG    ###################################################################################
###################################################################################################

def hdag_output_general(node_set, inp, taxId2seq, pars_weight="inf", bifurcate=False, adjust=False, true_tree_path=None):
    """
    Uses a generalized node support that can consider non-MP trees and weights them as a functions
    of their parsiomny score.

    Returns a list of tuples (clade_id, estimated_support, in_tree). The list is sorted by
    estimated support, and groups that have the same support are randomly shuffled.

    Args:
        node_set: Set of nodes (sets of taxa ids) that are in the true tree.
        inp: Otimized DAG or a protobuf file to the optimized DAG. WARNING: dag input may be
            mutated.
        taxId2Seq: Dicitonary of taxon IDs to the sequences they represent
        pars_weight: Constant to multiply by in the pscore_fn. If `inf`, then use uniform
            distribution over MP trees. Otherwise, support for edge (n1, n2) is computed as
            exp(-<pars_weight> * <pars(n1, n2)>)
        bifurcate: Whether to compute support with respect to all compatible bifurcating trees in
            the hDAG.
        adjust: Whether to use adjusted node support
    """

    print(f"=> Parsimony score weight: k={pars_weight}")

    larch_trim = False
    if larch_trim:
        executable = "/home/whowards/larch/larch/build/larch-usher"
        trimmed_path = ".".join(inp.split(".")[:-1]) + "_trimmed." + inp.split(".")[-1]
        print("Trimming:", inp, "-->", trimmed_path)
        args = [executable,
                "-i", f"{inp}",
                "-o", f"{trimmed_path}",
                "-l", "/fh/fast/matsen_e/whowards/hdag-benchmark/data/temp_logs",
                "--trim-max-parsimony",
                ]
        subprocess.run(args=args)
        inp = trimmed_path


    start = time.time()
    if isinstance(inp, str):
        dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(inp)
    else:
        dag = inp
    print(f"=> Loading DAG took {(time.time() - start)/60} minutes")

    start = time.time()
    # TODO: Complete and collapse in larch too
    # dag.make_complete()
    dag.convert_to_collapsed()
    print(f"=> Collapsing DAG took {(time.time() - start)/60} minutes")

    if isinstance(pars_weight, str):
        # This recovers uniform distribution over MP trees
        dag.trim_optimal_weight()
        pscore_fn = lambda n1, n2: 0
    else:
        pscore_fn = lambda n1, n2: -pars_weight * parsimony_utils.hamming_cg_edge_weight(n1, n2)

    num_trees = dag.count_trees()
    print("=> DAG contains", num_trees, "trees")

    if true_tree_path is not None:
        print("Getting true history...")
        true_history = get_history(true_tree_path, true_tree_path + ".fasta", dag)
        # Confirm that the two DAGs have the same leaf set
        assert {n.label for n in dag.get_leaves()} == {n.label for n in true_history.get_leaves()}

        rf_counts = dag.count_rf_distances(true_history)
        print("RF-distances to true tree:")
        print(rf_counts)
        # TODO: If the true tree is not contained in the DAG, you could add it and see how much
        #       support improves

    if bifurcate:
        def bifurcation_correction(node):
            if len(node.clades) > 2:
                return count_labeled_binary_topologies(len(node.clades))
            else:
                return 1
        dag.probability_annotate(lambda n1, n2: bifurcation_correction(n2), log_probabilities=False)
        log_prob = False
    else:
        # TODO: Refactor to use dag.uniform_annotate
        dag.probability_annotate(lambda n1, n2: pscore_fn(n1, n2), log_probabilities=True)
        log_prob = True

    if adjust:
        support = dag.adjusted_node_probabilities(log_probabilities=log_prob, collapse_key=lambda n: n.clade_union(), mut_type="site")
    else:
        support = dag.node_probabilities(log_probabilities=log_prob, collapse_key=lambda n: n.clade_union())

    seq2taxId = {v: k for k, v in taxId2seq.items()}

    return make_results_list(support, node_set, seq2taxId=seq2taxId, log_prob=log_prob)


@click.command("trim_thresholds")
@click.option('--tree_path', '-t', help='the newick file for the true tree.')
@click.option('--pb_file', '-i', help='the file containing inference results.')
@click.option('--output_path', '-o', help='the file to save output to.')
def trim_thresholds(tree_path, pb_file, output_path):
    """
    Given a path to a completed hDAG (e.g., data/A.2.2/1/results/historydag/final_opt_dag.pb)
    store results.pkls for various trimming strategies
    """

    fasta_path = tree_path + ".fasta"
    taxId2seq = hdag.parsimony.load_fasta(fasta_path)
    node_set = get_true_nodes(tree_path)

    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(pb_file)
    dag.make_complete()
    min_weight = dag.optimal_weight_annotate()

    strat_dict = {}

    strats = ["mp", 1.01, 1.05, 1.1, 1.2, "full"] # Proportion of parsimony to trim to
    for strat in strats:
        print(strat)
        if strats != "mp":
            dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(pb_file)
            dag.make_complete()
        
        if strat == "mp":
            dag.trim_optimal_weight()
        elif strat == "full":
            ... # Do no trimming
        else:
            max_weight = int(strat * min_weight)
            dag.trim_below_weight(max_weight=max_weight)
        
        trimmed_weights = dag.weight_count()
        counts = dag.count_nodes(collapse=True)
        total_trees = dag.count_trees()

        seq2taxId = {v: k for k, v in taxId2seq.items()}

        node2stats = {}
        counter = 0
        for node in counts:
            if len(node) == 0:  # UA node
                continue
            # Convert cg label to taxon id
            id_node = frozenset([seq2taxId[label.compact_genome.to_sequence()] for label in node])
            est_sup = counts[node] / total_trees
            node2stats[id_node] = (est_sup, id_node in node_set)
        # Get the support for all nodes in true tree
        for id_node in node_set:
            if id_node not in node2stats.keys():
                node2stats[id_node] = (0, True)
        stats_list =[(id_node, stats[0], stats[1]) for id_node, stats in node2stats.items()]
        random.shuffle(stats_list)
        stats_list.sort(key=lambda el: el[1])

        strat_dict[strat] = (trimmed_weights, stats_list)

    with open(output_path, "wb") as f:
        pickle.dump(strat_dict, f)


@click.command("test_pars_weights")
@click.option('--tree_path', '-t', help='the newick file for the true tree.')
@click.option('--pb_file', '-i', help='the file containing inference results.')
@click.option('--output_path', '-o', help='the file to save output to.')
def test_pars_weights(tree_path, pb_file, output_path):
    """
    Given a path to a completed hDAG (e.g., data/A.2.2/1/results/historydag/final_opt_dag.pb)
    store results.pkls for various parsiomny weightings
    """
    start = time.time()
    fasta_path = tree_path + ".fasta"
    taxId2seq = hdag.parsimony.load_fasta(fasta_path)
    node_set = get_true_nodes(tree_path)
    print("\n\nloading fasta and node_set took", time.time()-start, "seconds.\n\n")

    weight_dict = {}
    pweight = [-2, -1, -0.5, 0, 0.1, 0.5, 1, 2, "inf"] # Proportion of parsimony to trim to

    # print("Loading MAD...")
    # start = time.time()
    # dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(pb_file)
    # print(f"\t Took {(time.time() - start)/60} minutes")
    # TODO: ^ Can we save having to load this thing multiple times

    start = time.time()
    for p in pweight:
        stats_list = hdag_output_general(node_set, pb_file, taxId2seq, pars_weight=p, bifurcate=pweight == "bif")
        weight_dict[p] = (p, stats_list)
    print("\n\nEntire run of test_pars_weight took", time.time()-start, "seconds.\n")

    with open(output_path, "wb") as f:
        pickle.dump(weight_dict, f)


# TODO: Add a CLI arg that allows larch to collapse the DAG
@click.command("larch_usher")
@click.option('--executable', '-e', default='/home/whowards/larch/larch/build/larch-usher', help='path to pre-made larch-usher executable')
@click.option('--input', '-i', help='input tree or hdag. if tree, need refseqfile.')
@click.option('--refseqfile', '-r', default=None, help='number of .')
@click.option('--count', '-c', help='number of iterations.')
@click.option('--out_dir', '-o', help='the directory for where to store resulting dag protobufs.')
@click.option('--final_dag_name', '-f', default='final_opt_dag')
@click.option('--schedule', '-s', default="annealed")
@click.option('--log_dir', '-l')
@click.option('--pars_score', '-p', default=1)
@click.option('--node_score', '-n', default=1)
def larch_usher(executable, input, refseqfile, count, out_dir, final_dag_name, schedule, log_dir, pars_score, node_score):
    """Python CLI for driving larch-usher
    
    E.g.,
        larch_usher \
        -i /fh/fast/matsen_e/whowards/hdag-benchmark/data/A.2.5/1/results/historydag/opt_dag_3.pb \
        -c 2000 \
        -o /fh/fast/matsen_e/whowards/hdag-benchmark/data/A.2.5/1/results/historydag \
        -l /fh/fast/matsen_e/whowards/hdag-benchmark/data/A.2.5/1/results/historydag/opt_info
    """

    # E.g., results/historydag/
    os.chdir(f"{out_dir}")

    if int(count) <= 2:
        raise Exception("Not enough iterations")
        return


    # Cast a wide net by prioritizing new nodes only
    print("Running initial iterations of larch-usher...")
    subprocess.run(["mkdir", "-p", f"{log_dir}_1"])
    args = [executable,
            "-i", f"{input}",
            "-c", f"{round(int(count))}",
            "-o", f"{out_dir}/opt_dag_1.pb",
            "-l", f"{log_dir}_1",
            "--move-coeff-nodes", str(4),
            "--move-coeff-pscore", str(1),
            "--sample-any-tree"
            ]
    if refseqfile is not None:
        args.extend(["-r", refseqfile])
    subprocess.run(args=args)

    # Start considering parsimonious moves
    subprocess.run(["mkdir", "-p", f"{log_dir}_2"])
    args = [executable,
            "-i", f"{out_dir}/opt_dag_1.pb",
            "-c", f"{round(int(count))}",
            "-o", f"{out_dir}/opt_dag_2.pb",
            "-l", f"{log_dir}_2",
            "--move-coeff-nodes", str(1),
            "--move-coeff-pscore", str(1)
            ]
    subprocess.run(args=args)

    # Emphasize parsimony over new nodes
    subprocess.run(["mkdir", "-p", f"{log_dir}_3"])
    args = [executable,
            "-i", f"{out_dir}/opt_dag_2.pb",
            "-c", f"{round(int(count))}",
            "-o", f"{out_dir}/opt_dag_3.pb",
            "-l", f"{log_dir}_3",
            "--move-coeff-nodes", str(1),
            "--move-coeff-pscore", str(3),
            "--trim-max-parsimony"
            ]
    subprocess.run(args=args)

    print("Completing DAG...")
    start = time.time()
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(f"{out_dir}/opt_dag_3.pb")
    dag.convert_to_collapsed()
    dag.make_complete()
    dag.to_protobuf_file(f"{out_dir}/complete_opt_dag.pb")
    end = time.time()
    print(f"\tCompletion took {(end - start) / 60} minutes")

    subprocess.run(["mkdir", "-p", f"{log_dir}_complete"])
    args = [executable,
            "-i", f"{out_dir}/complete_opt_dag.pb",
            "-c", f"{round(int(count))}",
            "-o", f"{out_dir}/{final_dag_name}.pb",
            "-l", f"{log_dir}_complete",
            "--move-coeff-nodes", str(1),
            "--move-coeff-pscore", str(3),
            "--sample-any-tree",
            "--trim-max-parsimony"
            ]
    subprocess.run(args=args)


###################################################################################################
###   MrBayes   ###################################################################################
###################################################################################################


def get_mrbayes_trees(trees_file):
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
        for line in fh:
            if 'tree' in line.strip()[0:5]:
                nwk = line.strip().split(' ')[-1]
                tree = build_tree(nwk, translate_dict)
                for node in tree.iter_leaves():
                    node.name = translate_dict[node.name]
                # put original ambiguous sequences back on leaves
                yield tree

def mrbayes_output(node_set, tree_file, burnin=0.7, sample_freq=1000):
    """Same as hdag output, but for MrBayes"""
    for i, tree in enumerate(get_mrbayes_trees(tree_file)):
        continue
    num_trees = i
    num_burned = num_trees * burnin

    node2support = {}
    tree_count = 0
    for i, tree in enumerate(get_mrbayes_trees(tree_file)):
        if i * sample_freq < num_burned:
            continue
        tree_count += 1
        rerooted = reroot(tree.search_nodes(name="ancestral")[0])
        node2cu = {}
        curr_internal_name = 0
        for node in rerooted.traverse("postorder"):
            if node.is_leaf():
                cu = [node.name]
            else:
                node.name = f"internal{curr_internal_name}"
                curr_internal_name += 1
                cu = []
                for child in node.children:
                    cu.extend(list(node2cu[child.name]))
            
            cu = frozenset(cu)
            node2cu[node.name] = cu
            if cu not in node2support:
                node2support[cu] = 0
            # The unrooted topologies in trprobs are unique, so we don't need
            # to worry about uniqueness when rerooting them:
            node2support[cu] += 1

    for node, count in node2support.items():
        node2support[node] = count / tree_count

    return make_results_list(node2support, node_set)


###################################################################################################
###   UFBoot   ####################################################################################
###################################################################################################

def parse_splits(file_path):
    taxa_mapping = {}
    with open(file_path, "r") as f:
        in_taxa_mapping = False
        for line in f:
            if line.startswith(";"):
                break

            if in_taxa_mapping:
                pair = line.strip().split(" ")
                ufboot_id = pair[0][1:-1]
                taxa_id = pair[1][1:-1]
                taxa_mapping[ufboot_id] = taxa_id
            
            if line.startswith("TAXLABELS"):
                in_taxa_mapping = True
    
    with open(file_path, "r") as f:
        in_splits = False
        for line in f:
            if line.startswith(";"):
                in_splits = False

            if in_splits:
                pair = line.strip().split("\t")
                print(pair, "->", int(pair[0].strip()), pair[1].strip()[:-1].split(" "))
                num_trees = int(pair[0].strip())
                clade = [taxa_mapping[taxa] for taxa in pair[1].strip()[:-1].split(" ")]
                yield clade, num_trees

            if line.startswith("MATRIX"):
                in_splits = True

def ufboot_output(node_set, input_file):
    """Same as hdag output, but for UFBoot"""

    node2support = {}

    tree_count = 100
    for node, support in parse_splits(input_file):
        if support == 0:
            continue
        cu = frozenset(node)
        node2support[cu] = support / tree_count

    return make_results_list(node2support, node_set)


###################################################################################################
###    Hdag Diffused   ############################################################################
###################################################################################################
        
def nwk_sample_output(node_set, tree_file):
    """Same as hdag output, but for file of ete trees.
    
    Command to test this on A.2.2/1:
        python support_pipeline_scripts/cli.py save_supports \
        -m "random" \
        -t "/fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.132/1/simulation/collapsed_simulated_tree.nwk" \
        -i "/fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.132/1/results/random/random_trees.trees" \
        -o "/fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.132/1/results/random/results.pkl"
    """

    total_trees = 0
    node2count = {}
    with open(tree_file, "r") as f:
        for line in f:
            nw_str = line.strip()
            if len(nw_str) == 0:
                continue
            tree = ete.Tree(nw_str, format=9)

            node2cu = {}
            curr_internal_name = 0
            for node in tree.traverse("postorder"):
                if node.is_leaf():
                    cu = [node.name]
                else:
                    node.name = f"internal{curr_internal_name}"
                    curr_internal_name += 1
                    cu = []
                    for child in node.children:
                        cu.extend(list(node2cu[child.name]))
                
                cu = frozenset(cu)
                node2cu[node.name] = cu
                if cu not in node2count:
                    node2count[cu] = 0
                node2count[cu] += 1
            total_trees += 1

    node2support = {}
    for node, count in node2count.items():
        # NOTE: This is a debugging check
        if count > total_trees:
            print(f"=> Count is {count} with {total_trees} trees")
            print("Node")
            print(node)
        node2support[node] = count / total_trees

    return make_results_list(node2support, node_set)

@click.command("diffused_hdag_samples")
@click.option('--out_file', '-o', help='path to pre-made larch-usher executable')
@click.option('--input_dag', '-i', help='input MADAG. Does not trim or complete the given DAG..?')
@click.option('--num_samples', '-n', default=1000, help='number of ete trees to sample from the given DAG')
@click.option('--fasta_path', '-f', help='fasta path for given DAG that maps sequences to taxon id')
@click.option('--mut_rates_path', '-m', default=None, help='phastSim true mutation rates file')
@click.option('--branch_len', '-b', default=0.001, help='default branch length')
def diffused_hdag_samples(input_dag, out_file, num_samples, fasta_path, mut_rates_path, branch_len):
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(input_dag)
    dag.make_complete()
    dag.convert_to_collapsed()

    if mut_rates_path is not None:
        substitution_model = np.array(
                [
                    [-0.472, 0.039, 0.31, 0.123],
                    [0.14, -3.19, 0.022, 3.028],
                    [0.747, 0.113, -3.813, 2.953],
                    [0.056, 0.261, 0.036, -0.353]
                ]
            )
        site2rates = {}
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

        sampler = dag.diffused_tree_sampler(num_samples, load_fasta(fasta_path), lambda n: prob_muts(n, site2rates, branch_len=branch_len))
    else:
        sampler = dag.diffused_tree_sampler(num_samples, load_fasta(fasta_path), lambda ete_node: 0.02) 

    with open(out_file, "w") as fh:
        for _ in range(num_samples):
            ete_tree = next(sampler)
            # TODO: Figure out how to store mutations here?
            fh.write(f"{ete_tree.write(format=9)}\n")


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

def prob_muts_basic(ete_node, scale=0.01):
    muts = list(
                    cg_diff(
                        ete_node.up.compact_genome, ete_node.compact_genome
                    )
                )
    return scale * 1 / (len(muts)+1)




cli.add_command(larch_usher)
cli.add_command(save_supports)
cli.add_command(diffused_hdag_samples)


if __name__ == '__main__':
    cli()
