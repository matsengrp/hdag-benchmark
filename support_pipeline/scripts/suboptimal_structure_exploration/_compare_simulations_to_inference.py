import historydag as hdag
from predict_hotspots import predict_hotspots
from pathlib import Path
import json
from collections import Counter
import ete3
from io import BytesIO
import random
import subprocess
import os
random.seed(123)

# NOTE: Need to run the following in bash session
# export QT_QPA_PLATFORM=offscreen; export XDG_RUNTIME_DIR=/tmp/runtime-runner; export MPLBACKEND=agg; ml libGLU/9.0.1-GCCcore-10.2.0

from PyPDF2 import PdfReader, PdfWriter

outpath = Path("/fh/fast/matsen_e/whowards/hdag-benchmark/data/sub_struct")
outpath.mkdir(exist_ok=True)

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

def replace_site(site, newbase, sequence):
    idx = site - 1
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
            is_complex_pcm = diff_node.clade_union() == child_w_mut_cu

            if is_simple_pcm: # All children collapsed or doesn't count
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

path_prefix = "/fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models"
clade_names = [
    # "AY.34.2",
    # "AY.108",
    "AZ.3",
    "P.1.7",
    "A.2.5",
    "AY.74",
    "AY.87",
    "B.1.1.10"
]
hmut_params = [
    "_unrest_hmut_150_10",
    # "_unrest_hmut_200_10",
    # "gamma_rv_10",
    # "gamma_rv_20"
]
trials = range(1, 26)

param_trial_list = [(el1, el2, el3) for el1 in clade_names for el2 in hmut_params for el3 in trials]

# with open('support_pipeline/scripts/suboptimal_structure_exploration/to-consider.txt', 'w') as fh:
#     for clade_name, sim_type, trial in param_trial_list:
#         trial_path = f"{path_prefix}/{clade_name}/{sim_type}/{trial}"
#         fh.write(f"{trial_path}\n")

with open('support_pipeline/scripts/suboptimal_structure_exploration/to-consider.txt', 'r') as fh:
    inpaths = [Path(line.strip()) for line in fh if len(line) > 2 and line[0] != '#']


agg_size_data = []
agg_boundary_data = []
total_diffs = 0
duplicate_children = 0
for idx, inpath in enumerate(inpaths):
    idx += 42
    print(idx, inpath)
    dagpath = inpath / 'results/historydag/final_opt_dag_trimmed.pb'
    inpath = inpath / 'simulation'
    local_outpath = outpath / str(idx)
    local_outpath.mkdir(exist_ok=True)
    fasta = hdag.parsimony.load_fasta(inpath / 'ctree_with_refseq.fasta')
    variants_fasta = hdag.parsimony.load_fasta(inpath / 'collapsed_simulated_tree.nwk.variant_sites.fasta')
    etetree = ete3.Tree(str(inpath / 'collapsed_simulated_tree.nwk'), format=0)
    with open(inpath / 'collapsed_simulated_tree.nwk.variant_sites.txt', 'r') as fh:
        variant_sites = set(int(site) - 1 for site in fh.read().strip().split(' '))

    ancestral_seq = fasta["ancestral"]
    for node in etetree.traverse("preorder"):
        if node.is_root():
            parent_seq = ancestral_seq
        else:
            parent_seq = node.up.sequence
        if len(node.mutations) > 0:
            for mut in node.mutations.split("|"):
                mut = mut.strip()
                parent_seq = replace_site(int(mut[1:-1]), mut[-1], parent_seq)
        node.add_feature("sequence", parent_seq)

    # TODO: Fix this!
    #       /fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models/AY.87/_unrest_hmut_150_10/11
    #       ^-- This sim tree has a unifurcation
    try:
        sim_tree = hdag.history_dag_from_trees(
            [etetree],
            [],
            label_functions = {"sequence": lambda n: n.sequence},
            attr_func = lambda n: {"name": n.name if n.is_leaf() else "internal"},
        )
    except ValueError:
        continue

    try:
        opt_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(dagpath)
    except FileNotFoundError:
        continue

    # print("Loaded data, now completing and trimming")
    # opt_dag.convert_to_collapsed()
    # best_pars = opt_dag.optimal_weight_annotate()
    # print(f"trimmed to {opt_dag.count_trees()} trees, now relabeling")
    # remove invariant sites, convert to sequence DAG
    # new_cg_reference = subset_sequence(next(opt_dag.postorder()).label.compact_genome.reference, variant_sites)

    ref_seq = next(opt_dag.preorder(skip_ua_node=True)).label.compact_genome.reference
    sim_tree = hdag.mutation_annotated_dag.CGHistoryDag.from_history_dag(sim_tree, reference=ref_seq)
    sim_tree.convert_to_collapsed()
    print("relabeled, now adding sim tree")

    # TODO:  Can we also run sankoff
    # larch(sim_tree, str(local_outpath), count=1000)
    opt_sim_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(str(local_outpath) + "/trimmed_opt_dag.pb")

    opt_dag.merge([opt_sim_dag])
    opt_dag.make_complete()
    opt_dag.convert_to_collapsed()
    best_pars = opt_dag.optimal_weight_annotate()
    opt_dag.trim_optimal_weight()
    print(f"added sim tree. DAG has {opt_dag.count_trees()} trees")

    # TODO: Why does this fail for /fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models/AY.108/_unrest_hmut_150_10/18
    if set(sim_tree.get_leaves()) != set(opt_dag.get_leaves()):
        print(Warning(f"Skipping due to DAG and data not having same leaves!!! {dagpath}"))
        continue
    pars = sim_tree.optimal_weight_annotate()

    print(f"Best pars: {best_pars}\tSim pars: {pars}")



    # Add in leaf names in opt_dag, since they aren't preserved by protobuf:
    leaf_names = {n: n.attr["name"] for n in sim_tree.get_leaves()}
    for n in opt_dag.preorder():
        if n.is_leaf():
            n.attr = {"name": leaf_names[n]}
        else:
            n.attr = {"name": "internal"}

    hotspots, common_muts = predict_hotspots(opt_dag)
    opt_dag.trim_optimal_rf_distance(sim_tree)
    fitch_dag = opt_dag.copy()
    fitch_dag_pars = fitch_dag.optimal_weight_annotate()
    difference_filter = get_hamming_tree_compare(sim_tree)
    closest_pars_tree = fitch_dag[difference_filter][0]
    sim_tree.recompute_parents()
    closest_pars_tree.recompute_parents()
    #
    # We now have 'closest_pars_tree' and 'sim_tree' to compare, however we'd
    # like to.

    print("Finding suboptimal structures...")

    #############################################
    # def make_node(dag_node):
    #     """Given a DAG node, returns a tuple of clade union and CG"""
    #     return (dag_node.clade_union(), dag_node.label)
    # def get_up_edge(dag_node):
    #     """Given a DAG node, returns the edge to its parent."""
    #     if dag_node.is_root():
    #         return (None, make_node(dag_node))
    #     else:
    #         return (make_node(list(dag_node.parents)[0]), make_node(dag_node))
    # def get_cc(start_node, excluded_edges, visited_edges):
    #     # NOTE: Since we are starting from the root, we only need to consider children, not parents
    #     # visited = set()   # Unnecessary because graph is a directed tree
    #     cc = set()
    #     queue = [n for n in list(start_node.parents)[0].children() if get_up_edge(n) not in excluded_edges]
    #     while len(queue) > 0:
    #         node = queue[0]
    #         queue = queue[1:]
    #         edge = get_up_edge(node)
    #         cc.add(edge)
    #         visited_edges.add(edge)
    #         for child in node.children():
    #             child_edge = get_up_edge(child)
    #             if child_edge not in excluded_edges:
    #                 queue.append(child)
    #     return cc
    # pars_edges = set()
    # for node in closest_pars_tree.postorder():
    #     if not node.is_root():
    #         pars_edges.add(get_up_edge(node))
    # Find the size of all connected components using edges not in MP tree
    # visited_edges = set()
    # cc_counts = Counter()
    # ccs = set()
    # num_pars_edges = 0
    # for node in sim_tree.preorder():
    #     if node.is_root():
    #         continue
    #     edge = get_up_edge(node)
    #     if edge in pars_edges:
    #         num_pars_edges += 1
    #         continue
    #     elif edge in visited_edges:
    #         continue
    #     else:
    #         cc = get_cc(node, pars_edges, visited_edges)
    #         print(f"\tCC of size {len(cc)} starts at {len(list(node.parents)[0].clade_union())} sized node")
    #         # TODO: Get the height of this node
    #         cc_counts[len(cc)] += 1
    #         ccs.add(frozenset(cc))
    #     visited_edges.add(edge)
    # print("Distribution of connected components:")
    # print(cc_counts)
    # print("Number of pars edges:", num_pars_edges)
    # print("Number of edges by DAG:", sim_tree.num_edges())
    #############################################

        
    snode, pnode = get_difference_ancestor(sim_tree, closest_pars_tree)
    diffs, dups, pcm_clade_list = categorize_differences(closest_pars_tree, sim_tree)
    total_diffs += diffs
    duplicate_children += dups
    print(f"======\nTOTAL DIFFERENCES:\t{diffs}")
    print(f"TOTAL DUPLICATES:\t{dups}")
    print(f"\t--> Running tally: {duplicate_children} / {total_diffs}\n======\n")


    # Populate ete version of trees with relevant attributes
    same_nodes = set(sim_tree.preorder()) & set(closest_pars_tree.preorder())
    same_nodes = {node: idx for idx, node in enumerate(same_nodes)}
    kwargs = {"name_func": lambda n: n.attr.get("name", "internal"),
              "feature_funcs": {
                    "muts": lambda n: list(get_muts(n)),
                    "is_new": lambda n: n not in same_nodes,
                    "node_id": lambda n: str(same_nodes.get(n, "")),
                    # "in_cc": lambda n: [get_up_edge(n) in cc for cc in ccs]
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

print(f"Of {total_diffs} total topological differences, at least {duplicate_children} can be explained by duplicate child mutations")
# help from https://stackoverflow.com/a/74392556
# concatenate all plots into one pdf file.
result_document = outpath / 'new-comparison-plots.pdf'
writer = PdfWriter()
with open(result_document, 'wb') as fh:
    writer.write(fh)

# load examples into memory
examples = []
for filepath in list(outpath.glob('*/')):
    if filepath.is_dir() and int(str(filepath).split("/")[-1]) < len(inpaths):
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
    # height = simplotpage.mediabox.height + parsplotpage.mediabox.height
    # width = max(simplotpage.mediabox.width, parsplotpage.mediabox.width)
    # merged_page = PageObject.create_blank_page(None, width, height)
    # simplotpage.add_transformation(Transformation().scale(0.5).translate(0, 0))
    # merged_page.merge_page(simplotpage)
    # y = float(parsplotpage.mediabox.height)
    # parsplotpage.add_transformation(Transformation().scale(0.5).translate(0, y))
    # merged_page.merge_page(parsplotpage)
    writer.add_page(simplotpage)
    writer.add_page(parsplotpage)

with open(result_document, 'wb') as fh:
    writer.write(fh)


