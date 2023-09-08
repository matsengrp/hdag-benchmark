import click
import re
import ete3 as ete
import random
import pickle
import numpy as np
import os
import json

import random
import subprocess

# from newick_parser.tree_transformer import iter_nexus_trees
from historydag import parsimony_utils
from math import exp
import math
import time
from collections import Counter

import historydag as hdag
from historydag import parsimony_utils
from historydag.parsimony import parsimony_score, sankoff_upward, build_tree, sankoff_upward, load_fasta
from historydag.utils import count_labeled_binary_topologies

def get_true_nodes(tree_path, fasta=None):
    """ Given a newick file path, returns the set of all clades (as frozen sets of taxon ids)"""

    try:
        tree = ete.Tree(tree_path, format=3)
    except:
        print("\n --- Warning: Returning empty node set --- \n")
        return set()
    # NOTE: Removed this functionality for fast and non-silent failing
    # try:
    #     tree = ete.Tree(tree_path)
    # except:
    #     try:
    #         assert fasta is not None
    #         with open(tree_path, "r") as f:
    #             phastsim_tree = load_phastsim_newick(f.read())
    #             # TODO: Remove unifurcations and duplicate leaf nodes that have different names???
    #             tree = fix_phastsim_tree(phastsim_tree, fasta)
    #     except:
    #         print("\n --- Warning: Returning empty node set --- \n")
    #         return set()

    curr_internal_name = 0
    etenode2cu = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            cu = [node.name]
        else:
            # if (node.name) == 0:             # No name
            node.name = curr_internal_name
            curr_internal_name += 1
            cu = []
            for child in node.children:
                cu.extend(list(etenode2cu[child.name]))
        
        etenode2cu[node.name] = frozenset(cu)

    return set([v for k, v in etenode2cu.items()])


### Adapted from hdb collapse_tree.py #############################################################

# def load_phastsim_newick(newickstring):
#     """phastSim puts commas in node attributes, which ete can't parse"""
#     tree = ete.Tree()
#     current_node = tree

#     def internalnodecreate(sgen, part, current_node):
#         n = current_node.add_child()
#         return n

#     def attributecreate(sgen, part, current_node):
#         char = ''
#         current_node.name = part
#         while char != ']':
#             _, char = next(sgen)
#             part += char
#         strippedlist = part.split('{')[-1].split('}')[0]
#         if len(strippedlist) == 0:
#             mutset = set()
#         else:
#             mutset = set(strippedlist.split(','))
#         current_node.add_feature("mutations", mutset)
#         return current_node

#     def siblingcreate(sgen, part, current_node):
#         pnode = current_node.up
#         n = pnode.add_child()
#         return n

#     def moveup(sgen, part, current_node):
#         return current_node.up

#     def addbranchlength(sgen, part, current_node):
#         idx, char = next(sgen)
#         while newickstring[idx + 1] not in '([,):;':
#             idx, nchar = next(sgen)
#             char += nchar
#         current_node.dist = float(char)
#         return current_node

#     tokens = {'(': internalnodecreate,
#               '[': attributecreate,
#               ',': siblingcreate,
#               ')': moveup,
#               ':': addbranchlength,
#               ';': lambda _, __, cn: cn,}

#     sgen = enumerate(newickstring)
#     part = ''

#     while True:
#         try:
#             current_idx, char = next(sgen)
#             if char in tokens:
#                 current_node = tokens[char](sgen, part, current_node)
#                 part = ''
#             else:
#                 part += char
#         except StopIteration:
#             return tree
        
# def fix_phastsim_tree(etetree, fasta):
#     """Expects an ete tree with 'mutations' node attributes like that output by phastSim."""
#     def node_delete(node):
#         for child in node.children:
#             child.mutations.update(node.mutations)
#         node.delete(prevent_nondicotomic=False)

#     tree = etetree.copy()
#     # NOTE: Do everything except this
#     # # collapse zero length internal edges
#     # for node in tree.get_descendants():
#     #     if not node.is_leaf() and len(node.mutations) == 0:
#     #         node_delete(node)

#     # remove duplicate leaves which are grouped together
#     visited = set()
#     to_delete = []
#     for leaf in tree.get_leaves():
#         parent = leaf.up
#         if id(parent) not in visited:
#             visited.add(id(parent))
#             identical_children = [node for node in parent.children
#                                   if node.is_leaf() and len(node.mutations) == 0]
#             to_delete.extend(identical_children[1:])

#         # if leaf.name not in fasta.keys() and leaf.name not in to_delete:
#         #     to_delete.append(leaf)
    
#     for n in to_delete:
#         node_delete(n)

#     # remove unifurcation
#     for node in tree.get_descendants():
#         if len(node.children) == 1:
#             node_delete(node)

#     print(tree)

#     # check leaf sequences are unique
#     leaves = tree.get_leaves()
#     leaf_seqs = {fasta[node.name] for node in leaves}   # TODO: Why are there empty nodes in the tree?
#     n = len(leaves) - len(leaf_seqs)
#     if len(leaves) != len(leaf_seqs):
#         raise RuntimeError("Non-unique leaf sequences in modified tree")
#     return tree

###################################################################################################

# TODO: This could probably go back in the inference pipeline
def make_results_list(node2support, node_set, seq2taxId=None, log_prob=True):
    """Given a dictionary of node supports, and set of true nodes, returns a results list
    of triples, ordered by estimated support. `seq2taxId` and `log_prob` only need to be set
    if we are making a results list from the historydag inference"""

    print("Method supports", len(node2support), "nodes (with unique labels and clades)")

    if seq2taxId is not None:
        support = {}
        for node in node2support:
            if len(node) <= 1:  # UA node and leaves
                continue
            # Convert cg label to taxon id so that its compatible with true node ids
            id_node = frozenset([seq2taxId[label.compact_genome.to_sequence()] for label in node])
            est_sup = node2support[node]
            if log_prob:
                est_sup = exp(est_sup)
            support[id_node] = est_sup
        node2support = support

    # Construct results dict that maps nodes (frozen sets of taxon ids) to tuples of estimated
    # support and whether that node is in the true tree or not
    node2stats = {}

    # Get the support for all dag nodes
    for id_node, est_sup in node2support.items():
        if len(id_node) == 0:  # UA node
            continue
        
        # Clamp supports so that they are no larger than 1.0
        node2stats[id_node] = (min(est_sup, 1.0), id_node in node_set)
    
    # Get the support for all nodes in true tree
    for id_node in node_set:
        if id_node not in node2stats.keys():
            node2stats[id_node] = (0, True)

    print("Considering", len(node2stats), "(topological) nodes including those in true tree")
    stats_list =[(id_node, stats[0], stats[1]) for id_node, stats in node2stats.items()]
    random.shuffle(stats_list)
    stats_list.sort(key=lambda el: el[1])
    return stats_list


def reroot(new_root):
    """
    Edits the tree that the given node, new_root, is a part of so that it becomes the root.
    Returns pointer to the new root. Also, removes any unifurcations caused by edits.
    """
    node_path = [new_root]
    curr = new_root
    while not curr.is_root():
        node_path.append(curr.up)
        curr = curr.up

    root = node_path[-1]
    delete_root = len(root.children) <= 2
    
    while len(node_path) >= 2:
        curr_node = node_path[-1]
        curr_child = node_path[-2]
        curr_child.detach()
        curr_child.add_child(curr_node)
        node_path = node_path[:-1]
    if delete_root:
        root.delete()

    # Need to delete new root's child because it will be a unifurcation
    list(curr_child.children)[0].delete()

    return curr_child


def get_history(ete_tree_path, fasta_path, dag):
    """
    Given a path to a newick and DAG, returns a historyDAG containing a single history.
    """
    # TODO: Build ete tree with label attributes `sequence` containing the full nuc seqs at each node
    ete_tree = ete.Tree(ete_tree_path)
    id2seq = load_fasta(fasta_path)

    for node in ete_tree.traverse("postorder"):
        if node.is_leaf():
            node.add_feature("sequence", id2seq[node.name])
        else:
            prev_seq = None
            for child in node.children:
                new_seq = list(child.sequence)
                for mut in child.mutations.split("|"):
                    if len(mut) < 3:
                        continue
                    to = mut[0]
                    curr = mut[-1]
                    idx = int(mut[1:-1])-1
                    assert new_seq[idx] == curr
                    new_seq[idx] = to
                if prev_seq is not None:
                    assert new_seq == prev_seq
                prev_seq = new_seq
    
            node.add_feature("sequence", "".join(new_seq))

    seq_history = hdag.from_tree(ete_tree, label_features=["sequence"]) # ete_tree
    history = hdag.mutation_annotated_dag.CGHistoryDag.from_history_dag(seq_history, 
            reference=next(dag.postorder()).label.compact_genome.reference)

    return history


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


def larch(input, out_dir, count=1000, executable=None, move_coeff=1, pscore_coeff=3):
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
            "--move-coeff-nodes", str(move_coeff),
            "--move-coeff-pscore", str(pscore_coeff),
            "--sample-any-tree",
            ]
    with open(f"{log_dir}/mat_opt.log", "w") as f:
        subprocess.run(args=args, stdout=f, stderr=f)