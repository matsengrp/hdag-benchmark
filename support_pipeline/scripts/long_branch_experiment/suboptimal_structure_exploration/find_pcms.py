import historydag as hdag
from pathlib import Path
import ete3
from io import BytesIO
import random
import sys
import warnings
import subprocess
import os
import glob
import json
random.seed(123)



# NOTE: Need to run the following in bash session
# export QT_QPA_PLATFORM=offscreen; export XDG_RUNTIME_DIR=/tmp/runtime-runner; export MPLBACKEND=agg; ml libGLU/9.0.1-GCCcore-10.2.0


###################################################################################################
###   Helper methods                                                                            ###
###################################################################################################

def prune_tree(tree):
    """Remove children of the root whose subtrees contain no changes, leaving the root and (usually)
    a single child whose subtree contains all the changes"""
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node.modified_below = False
        else:
            node.modified_below = node.is_new or any(child.modified_below for child in node.children)

    to_remove = []
    for child in tree.children:
        if not child.modified_below:
            to_remove.append(child)
    for child in to_remove:
        tree.remove_child(child)

def larch(input, out_dir, count=1000, executable=None):
    """Python method for driving larch"""
    if executable is None:
        executable = 'larch-usher'

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
            "-o", f"{out_dir}/final_opt_dag_sim.pb",
            "-l", f"{log_dir}",
            "--move-coeff-nodes", str(1),
            "--move-coeff-pscore", str(5),
            ]
    with open(f"{log_dir}/sim_opt.log", "w") as f:
        subprocess.run(args=args, stdout=f, stderr=f)

            
def get_muts(node):
    if node.is_ua_node():
        yield from ()
    else:
        pnode = next(iter(node.parents))
        if pnode.is_ua_node():
            yield from ()
        else:
            for pbase, cbase, idx in hdag.compact_genome.cg_diff(pnode.label.compact_genome, node.label.compact_genome): #zip(pnode.label.sequence, node.label.sequence)):
                if pbase != cbase:
                    yield pbase + str(idx + 1) + cbase

def get_hamming_tree_compare(ref_tree):
    """Return addfuncdict implementing weight which is the sum of hamming distances
    between the sequences on a tree's nodes and the sequences on the corresponding nodes
    of in reference tree with the same topology. Assume all trees in DAG are equally near
    to the reference, topologically."""
    ref_dict = {n.clade_union(): n for n in ref_tree.preorder(skip_ua_node=True)}

    def edge_weight_func(n1, n2):
        return hdag.parsimony_utils.hamming_cg_edge_weight(n2, ref_dict.get(n2.clade_union(), n2))

    return hdag.utils.HistoryDagFilter(
        hdag.utils.AddFuncDict(
            {
                'start_func': lambda n: 0,
                'edge_weight_func': edge_weight_func,
                'accum_func': sum,
            },
            name='TreeSequenceDifference'
        ),
        min
    )

def replace_site(site, newbase, sequence, oldbase=None):
    idx = site - 1
    if oldbase is not None and sequence[idx] != oldbase:
        print(f"WARNING: previous base was not as expected (expected={oldbase} vs actual={sequence[idx]} at {site})")
    return sequence[:idx] + newbase + sequence[idx + 1:]

def subset_sequence(sequence, idxs):
    return ''.join(sequence[idx] for idx in sorted(idxs))

def categorize_differences(mp_tree, sim_tree, uncollapsed_sim_tree):
    """
    Finds all the differing clades between mp_tree and sim_tree. Returns the count of those clades
    and how many of those are due to parallel child mutations (PCM). Will also return a list of the
    clades (sets of leaf ids) that are due to PCMs.
    """
    mp_tree.recompute_parents()
    sim_tree.recompute_parents()
    uncollapsed_sim_tree.recompute_parents()
    mp_d = {n.clade_union(): n for n in mp_tree.preorder(skip_ua_node=True)}
    sim_d = {n.clade_union(): n for n in sim_tree.preorder(skip_ua_node=True)}
    uncollapsed_sim_d = {n.clade_union(): n for n in uncollapsed_sim_tree.preorder(skip_ua_node=True)}
    
    diff_clades = set(mp_d.keys()) - set(uncollapsed_sim_d.keys())

    total_diff_count = len(diff_clades)

    duplicate_child_mutation_count = 0
    pcm_clades = set()

    cgclade2idclade = {}
    for node in mp_tree.postorder():
        if node.is_leaf():
            idclade = frozenset([node.attr.get("name")])
        else:
            idclade = frozenset()
            for child in node.children():
                idclade = idclade.union(cgclade2idclade[child.clade_union()])
        cgclade2idclade[node.clade_union()] = idclade

    for diff_clade in diff_clades:
        # diff clade is the node in the MP tree which is not seen in the
        # simulated tree

        # case where multiple child mutations explain topological difference:
        diff_node = mp_d[diff_clade]
        parent_clade = next(iter(diff_node.parents)).clade_union()
        diff_muts = set(get_muts(diff_node))

        if parent_clade in sim_d:
            # Node in simulated history corresponding to the parent of the differing node
            sim_parent = sim_d[parent_clade]
            def check_child(child_clade):
                if child_clade not in sim_parent.clades:
                    return False
                else:
                    child_n = sim_parent.clades[child_clade].targets[0]
                    child_muts = set(get_muts(child_n))
                    if diff_muts.issubset(child_muts):
                        return True
                    else:
                        return False
            is_simple_pcm = all(check_child(clade) for clade in diff_node.clades)
            
            # Clade union of children in simulated tree with offending mutation
            child_w_mut_cu = set()
            for sim_child in sim_parent.children():
                child_muts = set(get_muts(sim_child))
                if diff_muts.issubset(child_muts) and \
                    len(sim_child.clade_union() & diff_node.clade_union()) > 0:
                        child_w_mut_cu = child_w_mut_cu.union(sim_child.clade_union())

            # If the resulting set of taxa is a subset of the node in the MP tree,
            # then we could reduce the number of mutaitons by making a new node
            # Therefore, it's a PCM.
            is_complex_pcm = diff_node.clade_union() == child_w_mut_cu

            if is_simple_pcm: # All children must be collapsed for it to count
                duplicate_child_mutation_count += 1
                pcm_clades.add(cgclade2idclade[diff_clade])
            elif is_complex_pcm:    # Detects different num children collapsed
                duplicate_child_mutation_count += 1
                pcm_clades.add(cgclade2idclade[diff_clade])
                print("\tCAUGHT NON-SIMPLE PCM!")

            # TODO: Can I remove this?
            # is_complex should capture simple PCMs too
            assert not (is_simple_pcm and not is_complex_pcm)


    return total_diff_count, duplicate_child_mutation_count, pcm_clades
                    
def get_difference_ancestor(t1, t2):
    """Get the most recent common ancestors of the subtrees containing all
    the differences between t1 and t2."""
    t1.recompute_parents()
    t2.recompute_parents()
    ns1 = set(t1.preorder())
    ns2 = set(t2.preorder())
    d1 = ns1 - ns2
    d2 = ns2 - ns1
    a1 = get_ancestor(d1, t1)
    a2 = get_ancestor(d2, t2)
    if a1 != a2:
        return (t1.dagroot, t2.dagroot)
    else:
        return (a1, a2)

def get_ancestor(node_set, history):
    """Get the MRCA-which-is-not-in-the-set of a set of nodes in a history"""
    clade = set()
    for node in node_set:
        clade.update(node.clade_union())
    for node in history.postorder():
        if node.clade_union().issuperset(clade):
            if node in node_set:
                return next(iter(node.parents))
            else:
                return node
    raise RuntimeError

###################################################################################################
###                                                                                             ###
###################################################################################################


path_prefix = "/n/fs/ragr-research/users/wh8114/prev/hdag-benchmark/data/sim_models"
outpath = Path("/n/fs/ragr-research/users/wh8114/prev/hdag-benchmark/data/sub_struct/results")
outpath.mkdir(exist_ok=True)

clade_name = sys.argv[1]
trial = sys.argv[2]
branch_multiplier = sys.argv[3]
inpath = Path(f"{path_prefix}/{clade_name}/{trial}/branch_multiplier_{branch_multiplier}")

print("Input path:", str(inpath))

# NOTE: Assumes that trim_dag.sh has been run on this clade/trial/branch-multiplier
results_path = inpath / 'results/historydag/'
file_paths = glob.glob(f"{results_path}/final_opt_dag_trimmed.pb.intermediate.*_dir")
assert len(file_paths) == 1
dagpath = f"{file_paths[0]}/final_opt_dag_trimmed.pb"
sim_inpath = inpath / 'simulation'
local_outpath = outpath / clade_name / str(trial) / f"branch_multiplier_{branch_multiplier}"
local_outpath.mkdir(exist_ok=True, parents=True)
fasta = hdag.parsimony.load_fasta(sim_inpath / 'ctree_with_refseq.fasta')
etetree = ete3.Tree(str(sim_inpath / 'collapsed_simulated_tree.nwk'), format=0)


# Reconstruct sequences on each node
ancestral_seq = fasta["ancestral"]
for node in etetree.traverse("preorder"):
    if node.is_root():
        assert len(node.mutations) == 0
        parent_seq = ancestral_seq
    else:
        parent_seq = node.up.sequence
    
    if len(node.mutations) > 0:
        mut_dict = {}
        mut_list = node.mutations.split("|")
        for mut_str in mut_list:
            parent_seq = replace_site(int(mut_str[1:-1]), mut_str[-1], parent_seq, oldbase=mut_str[0])
        
    node.add_feature("sequence", parent_seq)

seqs = set(fasta.values())

sim_history = hdag.history_dag_from_trees(
    [etetree],
    [],
    label_functions = {"sequence": lambda n: n.sequence},
    attr_func = lambda n: {"name": n.name if n.is_leaf() else "internal"},
)

opt_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(dagpath, compact_genomes=True)
opt_dag.trim_optimal_weight()
opt_dag = opt_dag.remove_label_fields(["node_id"])  # Makes the opt_dag comparable to sim_tree


for leaf in opt_dag.get_leaves():
    seq = leaf.label.compact_genome.to_sequence()
    if seq not in fasta.values():
        raise warnings.warn("DAG leaves are not the same as fasta leaves. Aborting.")

ref_seq = next(opt_dag.preorder(skip_ua_node=True)).label.compact_genome.reference
sim_history = hdag.mutation_annotated_dag.CGHistoryDag.from_history_dag(sim_history, reference=ref_seq)
sim_history_uncollapsed = sim_history.copy()
sim_history.convert_to_collapsed()

# sim_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(str(local_outpath) + "/sim_dag.pb", compact_genomes=True)
# sim_dag = sim_dag.remove_label_fields(["node_id"])  # Makes the opt_dag and sim_tree comparable
# sim_dag.trim_optimal_weight()
# sim_leaves = set(sim_dag.get_leaves())
# dag_leaves = set(opt_dag.get_leaves())
# for leaf_node in sim_leaves:
#     assert leaf_node.label.compact_genome.to_sequence() in fasta.values()

# NOTE: Once it gets to here, we know that the DAG and the tree are on the same set of sequences


print(f"Before adding sim tree, DAG trees: {opt_dag.count_trees()}")
opt_dag.merge([sim_history]) #[sim_dag])
print(f"After adding sim tree,  DAG trees: {opt_dag.count_trees()}")
# opt_dag.make_complete() # NOTE: This may be to slow
opt_dag.convert_to_collapsed()
best_pars = opt_dag.optimal_weight_annotate()
opt_dag.trim_optimal_weight()
print(f"After collapsing DAG trees:        {opt_dag.count_trees()}")

pars = sim_history.optimal_weight_annotate()

print(f"Best pars: {best_pars}\tSim pars: {pars}")


opt_dag.trim_optimal_rf_distance(sim_history)
fitch_dag = opt_dag.copy()
fitch_dag_pars = fitch_dag.optimal_weight_annotate()
difference_filter = get_hamming_tree_compare(sim_history)
closest_pars_tree = fitch_dag[difference_filter][0]
sim_history.recompute_parents()
closest_pars_tree.recompute_parents()

# We now have 'closest_pars_tree' and 'sim_history' to compare, however we'd like to.

print("Finding differing nodes...")
snode, pnode = get_difference_ancestor(sim_history, closest_pars_tree)
diffs, dups, pcm_clade_list = categorize_differences(closest_pars_tree, sim_history, sim_history_uncollapsed)
print(f"======")
print(f"TOTAL DIFFERENCES:\t{diffs}")
print(f"TOTAL DUPLICATES: \t{dups}")
print(f"======\n")


# TODO: Save total differences total duplicates and info
out_file_path = inpath / "results/historydag/pcm_info.json"
pcm_info_dict = {
    "diffs": diffs,
    "dups": dups
}
with open(out_file_path, 'w', encoding='utf-8') as f:
    json.dump(pcm_info_dict, f, ensure_ascii=False, indent=4)