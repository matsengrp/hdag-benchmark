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
from historydag.parsimony import parsimony_score, sankoff_upward, build_tree, sankoff_upward
from historydag.utils import count_labeled_binary_topologies

def get_true_nodes(tree_path):
    """ Given a newick file path, returns the set of all clades (as frozen sets of taxon ids)"""

    curr_internal_name = 0
    tree = ete.Tree(tree_path)
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

    print("DAG supports:", len(node2support))

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

        node2stats[id_node] = (est_sup, id_node in node_set)
    
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
