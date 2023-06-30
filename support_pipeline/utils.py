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

def get_true_nodes(tree_path):
    """ Given a newick file path, returns the set of all clades (as frozen sets of taxon ids)"""

    try:
        tree = ete.Tree(tree_path)
    except:
        print("\n --- Warning: Returning empty node set --- \n")
        return set()

    curr_internal_name = 0
    etenode2cu = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            cu = [node.name]
        else:
            if len(node.name) == 0:             # No name
                node.name = curr_internal_name
                curr_internal_name += 1
            cu = []
            for child in node.children:
                cu.extend(list(etenode2cu[child.name]))
        
        etenode2cu[node.name] = frozenset(cu)

    return set([v for k, v in etenode2cu.items()])

# TODO: This could probably go back in the inference pipeline
def make_results_list(node2support, node_set, seq2taxId=None, log_prob=True):
    """Given a dictionary of node supports, and set of true nodes, returns a results list
    of triples, ordered by estimated support. `seq2taxId` and `log_prob` only need to be set
    if we are making a results list from the historydag inference"""

    print("Method supports", len(node2support), "nodes")

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

    print("Considering", len(node2stats), "nodes")
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