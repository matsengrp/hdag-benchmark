import historydag as hdag
from predict_hotspots import predict_hotspots
from pathlib import Path
import ete3
from io import BytesIO
import random
import subprocess
import os
import sys
random.seed(123)

# NOTE: Need to run the following in bash session
# export QT_QPA_PLATFORM=offscreen; export XDG_RUNTIME_DIR=/tmp/runtime-runner; export MPLBACKEND=agg; ml libGLU/9.0.1-GCCcore-10.2.0

from PyPDF2 import PdfReader, PdfWriter

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
            "--move-coeff-nodes", str(1),
            "--move-coeff-pscore", str(3),
            ]
    with open(f"{log_dir}/mat_opt.log", "w") as f:
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

node_style_palette = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f']
def draw_tree_fragment_old(etetree, filename, title):
    """Draw a tree fragment, including leaf names, mutations along branches, with differing nodes' mutations
    highlighted in red. This assumes that the ete trees have attributes: is_new, node_id, muts, and in_cc.."""
    etetree = etetree.copy()
    for node in etetree.traverse("postorder"):
        nodestyle = ete3.NodeStyle()
        if node.is_leaf():
            node.first_leaf_name = int(node.name[1:])
        else:
            node.first_leaf_name = min(n.first_leaf_name for n in node.children)
        node.children.sort(key=lambda n: n.first_leaf_name)

        if any(node.in_cc):
            idx = ([i for i, val in enumerate(node.in_cc) if val][0] + 1) % len(node_style_palette)
            # nodestyle["size"] = 20
            nodestyle["bgcolor"] = node_style_palette[idx]
            node.set_style(nodestyle)
            for child in node.children:
                if any(child.in_cc):
                    continue
                else:
                    childstyle = ete3.NodeStyle()
                    childstyle["bgcolor"] = "white"
                    child.set_style(childstyle)


    ts = ete3.TreeStyle()
    ts.show_leaf_name = True
    ts.force_topology = True
    ts.title.add_face(ete3.TextFace(title, fsize=10), column=0)

    def layout_fn(node):
        if not node.is_leaf():
            ete3.faces.add_face_to_node(ete3.TextFace(node.node_id, fsize=8), node, column=1, position='branch-right')
        if not node.is_root():
            if any(node.in_cc):
                color = "Crimson"
            else:
                color = "Black"
            ete3.faces.add_face_to_node(ete3.TextFace(",".join(node.muts), fgcolor=color, fsize=8), node, column=1, position='branch-top')
        
    ts.rotation = 90
    ts.layout_fn = layout_fn

    # etetree.img_style["size"] = 30
    # etetree.img_style["fgcolor"] = "blue"

    etetree.render(filename, tree_style=ts)

def draw_tree_fragment(etetree, filename, title):
    """Draw a tree fragment, including leaf names, mutations along branches, with differing nodes' mutations
    highlighted in red."""
    etetree = etetree.copy()
    for node in etetree.traverse("postorder"):
        nodestyle = ete3.NodeStyle()
        if node.is_leaf():
            node.first_leaf_name = int(node.name[1:])
        else:
            node.first_leaf_name = min(n.first_leaf_name for n in node.children)
        node.children.sort(key=lambda n: n.first_leaf_name)

        if node.is_new:
            nodestyle["size"] = 10
            nodestyle["fgcolor"] = "Crimson"
            node.set_style(nodestyle)
            for child in node.children:
                if child.is_new or child.is_pcm:
                    continue
                else:
                    childstyle = ete3.NodeStyle()
                    # if child.node_id == "":     # child node clade is not in both mp and sim trees
                    #     child.node_id = random.randint(0, len(node_style_palette))
                    childstyle["bgcolor"] = node_style_palette[int(child.node_id) % len(node_style_palette)]
                    child.set_style(childstyle)
        if node.is_pcm:
            nodestyle["size"] = 20
            nodestyle["fgcolor"] = "DarkGreen"
            node.set_style(nodestyle)

    ts = ete3.TreeStyle()
    ts.show_leaf_name = True
    ts.force_topology = True
    ts.title.add_face(ete3.TextFace(title, fsize=10), column=0)

    def layout_fn(node):
        if not node.is_leaf():
            ete3.faces.add_face_to_node(ete3.TextFace(node.node_id, fsize=8), node, column=1, position='branch-right')
        if not node.is_root():
            if node.is_pcm:
                color = "DarkGreen"
            elif node.is_new:
                color = "Crimson"
            else:
                color = "Black"
            ete3.faces.add_face_to_node(ete3.TextFace(",".join(node.muts), fgcolor=color, fsize=8), node, column=1, position='branch-top')

    ts.rotation = 90
    ts.layout_fn = layout_fn

    # etetree.img_style["size"] = 30
    # etetree.img_style["fgcolor"] = "blue"

    etetree.render(filename, tree_style=ts)

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

def categorize_differences(mp_tree, sim_tree):
    """
    Finds all the differing clades between mp_tree and sim_tree. Returns the count of those clades
    and how many of those are due to parallel child mutations (PCM). Will also return a list of the
    clades (sets of leaf ids) that are due to PCMs.
    """
    mp_tree.recompute_parents()
    sim_tree.recompute_parents()
    mp_d = {n.clade_union(): n for n in mp_tree.preorder(skip_ua_node=True)}
    sim_d = {n.clade_union(): n for n in sim_tree.preorder(skip_ua_node=True)}
    diff_clades = set(mp_d.keys()) - set(sim_d.keys())

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

            # DEBUG stuff
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


path_prefix = "/fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models"
outpath = Path("/fh/fast/matsen_e/whowards/hdag-benchmark/data/sub_struct")
outpath.mkdir(exist_ok=True)
clade_name = sys.argv[1]
trial = sys.argv[2]

sim_type = "gamma_10_hmut_50"
inpath = Path(f"{path_prefix}/{clade_name}/{sim_type}/{trial}")

print("Input path:", str(inpath))

dagpath = inpath / 'results/historydag/final_opt_dag_trimmed.pb'
inpath = inpath / 'simulation'
local_outpath = outpath / clade_name / str(trial)
local_outpath.mkdir(exist_ok=True)
fasta = hdag.parsimony.load_fasta(inpath / 'ctree_with_refseq.fasta')
etetree = ete3.Tree(str(inpath / 'collapsed_simulated_tree.nwk'), format=0)


# Reconstruct sequences on each node
ancestral_seq = fasta["ancestral"]  # TODO: Try refseq.fasta sequence
for node in etetree.traverse("preorder"):
    if node.is_root():
        assert len(node.mutations) == 0
        parent_seq = ancestral_seq
    else:
        parent_seq = node.up.sequence
    
    if len(node.mutations) > 0:
        # NOTE: Potentially need to resolve multiple mutations at a single site
        # - Currently just trying applying the mutations in order
        mut_dict = {}
        mut_list = node.mutations.split("|")
        for mut_str in mut_list:
            parent_seq = replace_site(int(mut_str[1:-1]), mut_str[-1], parent_seq, oldbase=mut_str[0])
        
    node.add_feature("sequence", parent_seq)

seqs = set(fasta.values())
# NOTE: Debug check. can commetn out or delete
for leaf in etetree.get_leaves():
    if leaf.sequence not in seqs:
        print("Mutations on naughty leaf:")
        mut_dict = {}
        curr = leaf
        while curr.up is not None:
            print("\t", curr.mutations.strip().split("|"))
            if len(curr.mutations) > 0:
                for mut in curr.mutations.strip().split("|"):
                    idx = int(mut[1:-1])
                    from_nuc = mut[0]
                    to_nuc = mut[-1]
                    if idx not in mut_dict:
                        mut_dict[idx] = []
                    mut_dict[idx].append((from_nuc, to_nuc))
            curr = curr.up
        
        for idx, pair_list in mut_dict.items():
            print(f"{idx}\t{pair_list}")

sim_history = hdag.history_dag_from_trees(
    [etetree],
    [],
    label_functions = {"sequence": lambda n: n.sequence},
    attr_func = lambda n: {"name": n.name if n.is_leaf() else "internal"},
)

opt_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(dagpath)

for leaf in opt_dag.get_leaves():
    seq = leaf.label.compact_genome.to_sequence()
    if seq not in fasta.values():
        raise Warning("DAG leaves are not the same as fasta leaves. Aborting.")

ref_seq = next(opt_dag.preorder(skip_ua_node=True)).label.compact_genome.reference
sim_history = hdag.mutation_annotated_dag.CGHistoryDag.from_history_dag(sim_history, reference=ref_seq)
sim_history.convert_to_collapsed()
print("relabeled, now adding sim tree")


sim_leaves = set(sim_history.get_leaves())
dag_leaves = set(opt_dag.get_leaves())

if sim_leaves != dag_leaves:

    ###################################################################################################
    ##### DEGUB DOJO
    ###################################################################################################

    sim_ref = next(sim_history.preorder(skip_ua_node=True)).label.compact_genome.reference
    opt_ref = next(opt_dag.preorder(skip_ua_node=True)).label.compact_genome.reference
    if sim_ref != opt_ref:
        print("sim:", sim_ref, "\nopt:", opt_ref)
    else:
        print("reference sequences are equal")
    
    sim_diff = sim_leaves.difference(dag_leaves)
    dag_diff = dag_leaves.difference(sim_leaves)

    print("Number of differing leaves is", len(sim_diff), len(dag_diff))
    
    sim_mut_sets = [set([(idx, mut) for idx, mut in leaf.label.compact_genome.mutations.items()]) for leaf in sim_diff]
    dag_mut_sets = [set([(idx, mut) for idx, mut in leaf.label.compact_genome.mutations.items()]) for leaf in dag_diff]
    sim_mut_sets.sort(key = lambda n: len(n))
    dag_mut_sets.sort(key = lambda n: len(n))
    for i, (sim_mut_set, dag_mut_set) in enumerate(zip(sim_mut_sets, dag_mut_sets)):
        sim_mut_diff = sim_mut_set.difference(dag_mut_set)
        dag_mut_diff = dag_mut_set.difference(sim_mut_set)

        print(i)
        print(sim_mut_set)
        print(dag_mut_set)
        print("\tsim mutation diffs:", sim_mut_diff)
        print("\tdag mutation diffs:", dag_mut_diff)

    # NOTE: The DAG leaves are always correct.
    print("Sim leaves in fasta:")
    for i, cg in enumerate(sim_diff):
        seq = cg.label.compact_genome.to_sequence()
        if seq in fasta.values():
            print(i, cg.label.compact_genome)
    print("DAG leaves in fasta:")
    for i, cg in enumerate(dag_diff):
        seq = cg.label.compact_genome.to_sequence()
        if seq in fasta.values():
            print(i, cg.label.compact_genome)
    
    print("----------")
    print()
    ###################################################################################################
    ###################################################################################################

    raise Warning(f"Skipping due to DAG and data not having same leaves!!! {dagpath}")

larch(sim_history, str(local_outpath), count=1000)
opt_sim_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(str(local_outpath) + "/trimmed_opt_dag.pb")

opt_dag.merge([opt_sim_dag])
opt_dag.make_complete()
opt_dag.convert_to_collapsed()
best_pars = opt_dag.optimal_weight_annotate()
opt_dag.trim_optimal_weight()
print(f"added sim tree. DAG has {opt_dag.count_trees()} trees")

pars = sim_history.optimal_weight_annotate()

print(f"Best pars: {best_pars}\tSim pars: {pars}")


# Add in leaf names in opt_dag, since they aren't preserved by protobuf:
leaf_names = {n: n.attr["name"] for n in sim_history.get_leaves()}
for n in opt_dag.preorder():
    if n.is_leaf():
        n.attr = {"name": leaf_names[n]}
    else:
        n.attr = {"name": "internal"}

hotspots, common_muts = predict_hotspots(opt_dag)
opt_dag.trim_optimal_rf_distance(sim_history)
fitch_dag = opt_dag.copy()
fitch_dag_pars = fitch_dag.optimal_weight_annotate()
difference_filter = get_hamming_tree_compare(sim_history)
closest_pars_tree = fitch_dag[difference_filter][0]
sim_history.recompute_parents()
closest_pars_tree.recompute_parents()
#
# We now have 'closest_pars_tree' and 'sim_history' to compare, however we'd
# like to.

print("Finding differing nodes...")
snode, pnode = get_difference_ancestor(sim_history, closest_pars_tree)
diffs, dups, pcm_clade_list = categorize_differences(closest_pars_tree, sim_history)
print(f"======\nTOTAL DIFFERENCES:\t{diffs}")
print(f"TOTAL DUPLICATES:\t{dups}")
print(f"======\n")

# Populate ete version of trees with relevant attributes
same_nodes = set(sim_history.preorder()) & set(closest_pars_tree.preorder())
same_nodes = {node: idx for idx, node in enumerate(same_nodes)}
kwargs = {"name_func": lambda n: n.attr.get("name", "internal"),
            "feature_funcs": {
                "muts": lambda n: list(get_muts(n)),
                "is_new": lambda n: n not in same_nodes,
                "node_id": lambda n: str(same_nodes.get(n, "")),
            }
            }
ssubtree = snode.to_ete_recursive(**kwargs)
psubtree = pnode.to_ete_recursive(**kwargs)

node2clade = {}
for tree in [psubtree, ssubtree]:
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            cu = frozenset([node.name])
        else:
            cu = frozenset()
            for child in node.children:
                cu = cu.union(node2clade[child])
        node2clade[node] = cu
        node.add_feature("is_pcm", node2clade[node] in pcm_clade_list)


copylist = [ssubtree.copy(), psubtree.copy()]
for ctree in copylist:
    prune_tree(ctree)
n_nodes = min(len(list(ctree.traverse())) for ctree in copylist)
if n_nodes > 2:
    ssubtree, psubtree = copylist
node_ids = {}
draw_tree_fragment(ssubtree, str(local_outpath / "sim_tree.pdf"), "simulated " + '/'.join(str(inpath).split('/')[-6:-1]) + f" pscore {pars}\nhotspots: {hotspots.most_common(10)}\nmost common mutations: {common_muts.most_common(10)}\ndiffs: {diffs} dups: {dups}")
draw_tree_fragment(psubtree, str(local_outpath / "pars_tree.pdf"), "optimal " + '/'.join(str(inpath).split('/')[-6:-1]) + f" pscore {best_pars}")

# help from https://stackoverflow.com/a/74392556
# concatenate all plots into one pdf file.
result_document = local_outpath / 'new-comparison-plots.pdf'
writer = PdfWriter()
with open(result_document, 'wb') as fh:
    writer.write(fh)

# load examples into memory
examples = []
for filepath in list(local_outpath.glob('*/')):
    if filepath.is_dir():
        try:
            with (open(filepath / "sim_tree.pdf", 'rb') as simplot,
                  open(filepath / "pars_tree.pdf", 'rb') as parsplot):
                examples.append((BytesIO(simplot.read()), BytesIO(parsplot.read()), str(filepath)))
        except FileNotFoundError:
            continue

for simplot, parsplot, filepath in examples:
    simplotreader = PdfReader(simplot)
    simplotpage = simplotreader.pages[0]
    parsplotreader = PdfReader(parsplot)
    parsplotpage = parsplotreader.pages[0]
    writer.add_page(simplotpage)
    writer.add_page(parsplotpage)

with open(result_document, 'wb') as fh:
    writer.write(fh)


