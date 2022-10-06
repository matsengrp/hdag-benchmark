"""Load tree with simulated sequences and collapse"""

def collapse_tree(etetree, fasta):
    """Expects an ete tree with 'mutations' node attributes like that output by phastSim."""
    tree = etetree.copy()
    # collapse zero length internal edges
    for node in tree.get_descendants():
        if not node.is_leaf() and len(node.mutations) == 0:
            node.delete(prevent_nondicotomic=False)

    # remove duplicate leaves which are grouped together
    visited = set()
    for leaf in tree.get_leaves():
        parent = leaf.up
        if id(parent) not in visited:
            visited.add(id(parent))
            identical_children = [node for node in parent.children
                                  if node.is_leaf() and len(node.mutations) == 0]
            for n in identical_children[1:]:
                n.delete(prevent_nondicotomic=False)

    # remove unifurcation
    for node in tree.get_descendants():
        if len(node.children) == 1:
            node.delete(prevent_nondicotomic=False)

    # check leaf sequences are unique
    leaves = tree.get_leaves()
    leaf_seqs = {fasta[node.name] for node in leaves}
    if len(leaves) != len(leaf_seqs):
        print("WARNING: non-unique leaf sequences in modified tree")

    return tree
