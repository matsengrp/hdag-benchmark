import gctree
import historydag as hdag
from pathlib import Path
from typing import NamedTuple

def aggregate_dnapars_trees(outfiles, outgroup, abundance_file):
    dag = None
    for outfile in outfiles:
        trees = gctree.phylip_parse.parse_outfile(outfile, root=outgroup, abundance_file=abundance_file)
        partial_dag = gctree.branching_processes._make_dag(trees, from_copy=False)
        for node in partial_dag.preorder(skip_root=True):
            if not node.is_leaf() and node.label.abundance > 0:
                node.label = type(node.label)(node.label.sequence, 0)
        if dag is None:
            dag = partial_dag
        else:
            dag.merge(partial_dag)
    return dag

def merge_dags(dags):
    dag = None
    for partial_dag in dags:
        print("Parsimony score of best tree in this dag is ", partial_dag.optimal_weight_annotate())
        if dag is None:
            dag = partial_dag
        else:
            dag.merge(partial_dag)

    print("final DAG stats:")
    dag.summary()
    return dag
