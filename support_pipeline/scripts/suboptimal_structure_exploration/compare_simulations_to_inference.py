import historydag as hdag
from predict_hotspots import predict_hotspots
from pathlib import Path
import json
import ete3
from ete3 import TreeStyle
from io import BytesIO

from PyPDF2 import PdfReader, PdfWriter, PageObject, Transformation

outpath = Path("new-suboptimal_topologies")
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
            
def get_muts(node):
    if node.is_ua_node():
        yield from ()
    else:
        pnode = next(iter(node.parents))
        if pnode.is_ua_node():
            yield from ()
        else:
            for idx, (pbase, cbase) in enumerate(zip(pnode.label.sequence, node.label.sequence)):
                if pbase != cbase:
                    yield pbase + str(idx + 1) + cbase

node_style_palette = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f']
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
                if child.is_new:
                    continue
                else:
                    childstyle = ete3.NodeStyle()
                    childstyle["bgcolor"] = node_style_palette[int(child.node_id) % len(node_style_palette)]
                    child.set_style(childstyle)



    ts = ete3.TreeStyle()
    ts.show_leaf_name = True
    ts.force_topology = True
    ts.title.add_face(ete3.TextFace(title, fsize=10), column=0)

    def layout_fn(node):
        if not node.is_leaf():
            ete3.faces.add_face_to_node(ete3.TextFace(node.node_id, fsize=8), node, column=1, position='branch-right')
        if not node.is_root():
            if node.is_new:
                color = "Crimson"
            else:
                color = "Black"
            ete3.faces.add_face_to_node(ete3.TextFace(",".join(node.muts), fgcolor=color, fsize=8), node, column=1, position='branch-top')

    ts.rotation = 90
    ts.layout_fn = layout_fn
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
    mp_tree.recompute_parents()
    sim_tree.recompute_parents()
    mp_d = {n.clade_union(): n for n in mp_tree.preorder(skip_ua_node=True)}
    sim_d = {n.clade_union(): n for n in sim_tree.preorder(skip_ua_node=True)}
    diff_clades = set(mp_d.keys()) - set(sim_d.keys())

    total_diff_count = len(diff_clades)

    duplicate_child_mutation_count = 0

    for diff_clade in diff_clades:
        # diff clade is the node in the MP tree which is not seen in the
        # simulated tree

        # case where multiple child mutations explain topological difference:
        diff_node = mp_d[diff_clade]
        parent_clade = next(iter(diff_node.parents)).clade_union()
        diff_muts = set(get_muts(diff_node))
        if parent_clade in sim_d:
            sim_parent = sim_d[parent_clade]
            # TODO: check if mutations on parent edge of diff_clade node are
            # duplicated on all edges pointing to nodes with child clades of
            # diff_clade node.
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
            if all(check_child(clade) for clade in diff_node.clades):
                duplicate_child_mutation_count += 1
    print(f"{total_diff_count} total differences, {duplicate_child_mutation_count} are duplicate child mutations")
    return total_diff_count, duplicate_child_mutation_count
                    



def get_difference_ancestor(t1, t2):
    """Get the latest common ancestors of the subtrees containing all
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


# Get list of paths for relevant simulations
with open('to-consider3.txt', 'r') as fh:
    inpaths = [Path(line.strip()) for line in fh]


agg_size_data = []
agg_boundary_data = []
total_diffs = 0
duplicate_children = 0
for idx, inpath in enumerate(inpaths[:20]):
    print(idx, inpath)
    dagpath = inpath / 'results/historydag/final_opt_dag_trimmed.pb'
    inpath = inpath / 'simulation'
    local_outpath = outpath / str(idx)
    local_outpath.mkdir(exist_ok=True)
    # # Get tree stats for reference
    # with open(inpath / 'tree_stats.json', 'r') as fh:
    #     stats = json.loads(fh.read())
    # # pars, best_pars = stats["pars_score"], stats["max_score_data"]
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


    sim_tree = hdag.history_dag_from_trees(
        [etetree],
        [],
        label_functions = {"sequence": lambda n: n.sequence},
        attr_func = lambda n: {"name": n.name if n.is_leaf() else "internal"},
    )



    try:
        opt_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(dagpath)
    except FileNotFoundError:
        continue
    print("Loaded data, now completing and trimming")
    opt_dag.make_complete()
    copy_opt_dag = opt_dag.copy()
    opt_dag.trim_optimal_weight()
    opt_dag.convert_to_collapsed()
    best_pars = opt_dag.optimal_weight_annotate()
    print("trimmed, now relabeling")
    # remove invariant sites, convert to sequence DAG
    new_cg_reference = subset_sequence(next(opt_dag.postorder()).label.compact_genome.reference, variant_sites)

    
    ref_seq = next(opt_dag.preorder(skip_ua_node=True)).label.compact_genome.reference

    sim_tree = hdag.mutation_annotated_dag.CGHistoryDag.from_history_dag(sim_tree, reference=ref_seq)
    assert set(sim_tree.get_leaves()) == set(opt_dag.get_leaves())
    pars = sim_tree.optimal_weight_annotate()

    # def sim_tree_relabel_func(node):
    #     it = node.label.compact_genome.subset_sites(variant_sites, one_based=False, new_reference=new_cg_reference)
    #     cgdict = {}
    #     for idx, (variantsite, base) in enumerate(zip(variant_sites, node.label.sequence))
    #     cg = hdag.compact_genome.CompactGenome({variant_sites[idx]: (ref_seq[variant_sites[idx] - 1], )})
    #     return [it]

    # print("relabel 2... now adding names")
    # assert False
    # opt_dag = opt_dag.update_label_fields(["compact_genome"], relabel_func)
    # print("relabel 3... now adding names")
    # opt_dag = hdag.sequence_dag.SequenceHistoryDag.from_history_dag(opt_dag)
    # print("relabeled... now adding names")
    # assert len(list((opt_dag | sim_tree).get_leaves())) ==  len(list(opt_dag.get_leaves()))

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
    ###### NOTICE I commented this out, which may cause problems?
    # hdag.parsimony.sankoff_downward(fitch_dag, partial_node_list = list(fitch_dag.postorder())[:-2])
    # assert sim_tree.optimal_weight_annotate() == pars
    fitch_dag_pars = fitch_dag.optimal_weight_annotate()
    # if fitch_dag_pars > best_pars:
    #     print(str(inpath), f"larch did not find as good a parsimony score as dnapars: {fitch_dag_pars} > {best_pars}")
    # elif best_pars > fitch_dag_pars:
    #     print(str(inpath), f"larch found better parsimony score than dnapars: {fitch_dag_pars} < {best_pars}")
    # best_pars = fitch_dag_pars


    difference_filter = get_hamming_tree_compare(sim_tree)
    closest_pars_tree = fitch_dag[difference_filter][0]
    #
    # We now have 'closest_pars_tree' and 'sim_tree' to compare, however we'd
    # like to.


    sim_tree.recompute_parents()
    closest_pars_tree.recompute_parents()


    # TODO: Find all subotpimal structures
    # - Get set of edges in closest_tree (ie, pair of parent clade + seq and child clade + seq)
    # - Find all connected components that consist of edges in this set

    def get_key(node):
        return (node.clade_union(), node.label)

    def get_node_dict(tree):
        return {get_key(n): n for n in tree.preorder(skip_ua_node=True)}


    sim_nodes = get_node_dict(sim_tree)
    infer_nodes = get_node_dict(closest_pars_tree)
    
    def examine_connected_difference_components(t1_nodedict, t2_nodedict):
        t1_anomalous = set(t1_nodedict.keys()) - set(t2_nodedict.keys())
        t2_anomalous = set(t2_nodedict.keys()) - set(t1_nodedict.keys())
        print(len(t1_anomalous), len(t2_anomalous))

        def traverse_outward(node, last_visited=None, leaf_condition = lambda n: False):
            if last_visited is None:
                print("starting traversal")
            assert len(node.parents) <= 1

            yield (node, leaf_condition(node))
            # print(get_key(node))
            if leaf_condition(node):
                print(True)
            if not leaf_condition(node):
                if len(node.parents) > 0:
                    pnode = next(iter(node.parents))
                    if last_visited is None or pnode != last_visited:
                        yield from traverse_outward(pnode, last_visited=node)

                for cnode in node.children():
                    if last_visited is None or cnode != last_visited:
                        yield from traverse_outward(pnode, last_visited=node)

        print("building chunks")
        chunks1 = []
        for node_key in t1_anomalous:
            traversal = list(traverse_outward(t1_nodedict[node_key], leaf_condition = lambda n: get_key(n) not in t1_anomalous))
            chunks1.append(([get_key(tnode) for tnode, leaf in traversal if leaf],
                           [get_key(tnode) for tnode, leaf in traversal if not leaf]))

        chunks2 = []
        for node_key in t2_anomalous:
            traversal = list(traverse_outward(t2_nodedict[node_key], leaf_condition = lambda n: get_key(n) not in t2_anomalous))
            chunks2.append(({get_key(tnode) for tnode, leaf in traversal if leaf},
                           {get_key(tnode) for tnode, leaf in traversal if not leaf}))

        # Now match by inclusion of boundary nodes, using clade union and sequence
        unique_chunks = dict()
        opaired = set()

        print("matching chunks, sim tree has", len(chunks1), "and optimal tree has", len(chunks2), "chunks")
        for boundary, interior in chunks1:
            paired = False
            for oboundary, ointerior in chunks2:
                if boundary.issubset(oboundary):
                    unique_chunks[oboundary] = (interior, ointerior)
                    opaired.add(frozenset(oboundary))
                elif oboudary.issubset(boundary):
                    unique_chunks[boundary] = (interior, ointerior)
                    opaired.add(frozenset(oboundary))
                else:
                    continue
                paired = True
                break
            if not paired:
                unique_chunks[boundary] = (interior, set())
        for boundary, interior in chunks2:
            if frozenset(boundary) not in opaired:
                unique_chunks[boundary] = (set(), interior)

        size_data = [len(interior) + len(ointerior) for interior, ointerior in unique_chunks.values()]
        boundary_data = [len(key) for key in unique_chunks.keys()]
        print("tree has", len(unique_chunks), "suboptimal substructures with boundary sizes", boundary_data, "and sizes", size_data)
        agg_size_data.append(size_data)
        agg_boundary_data.append(boundary_data)

    "starting chunk building!"
    examine_connected_difference_components(sim_nodes, infer_nodes)
        
    # snode, pnode = get_difference_ancestor(sim_tree, closest_pars_tree)
    # diffs, dups = categorize_differences(closest_pars_tree, sim_tree)
    # total_diffs += diffs
    # duplicate_children += dups
    # same_nodes = set(sim_tree.preorder()) & set(closest_pars_tree.preorder())
    # same_nodes = {node: idx for idx, node in enumerate(same_nodes)}
    # kwargs = {"name_func": lambda n: n.attr.get("name", "internal"),
        #       "feature_funcs": {"sequence": lambda n: None if n.is_ua_node() else n.label.sequence,
        #                         "muts": lambda n: list(get_muts(n)),
        #                         "is_new": lambda n: n not in same_nodes,
        #                         "node_id": lambda n: str(same_nodes.get(n, ""))},
        #       }
    # ssubtree = snode.to_ete_recursive(**kwargs)
    # psubtree = pnode.to_ete_recursive(**kwargs)
    # copylist = [ssubtree.copy(), psubtree.copy()]
    # for ctree in copylist:
        # prune_tree(ctree)
    # n_nodes = min(len(list(ctree.traverse())) for ctree in copylist)
    # if n_nodes > 2:
        # ssubtree, psubtree = copylist
    # node_ids = {}
    # # draw_tree_fragment(ssubtree, str(local_outpath / "sim_tree.pdf"), "simulated " + '/'.join(str(inpath).split('/')[-3:-1]))
    # # draw_tree_fragment(psubtree, str(local_outpath / "pars_tree.pdf"), "optimal " + '/'.join(str(inpath).split('/')[-3:-1]))
    # draw_tree_fragment(ssubtree, str(local_outpath / "sim_tree.pdf"), "simulated " + '/'.join(str(inpath).split('/')[-3:-1]) + f" pscore {pars}\nhotspots: {hotspots.most_common(10)}\nmost common mutations: {common_muts.most_common(10)}")
    # draw_tree_fragment(psubtree, str(local_outpath / "pars_tree.pdf"), "optimal " + '/'.join(str(inpath).split('/')[-3:-1]) + f" pscore {best_pars}")

print(f"Of {total_diffs} total topological differences, at least {duplicate_children} can be explained by duplicate child mutations")
# help from https://stackoverflow.com/a/74392556
# concatenate all plots into one pdf file.
result_document = outpath / 'new-comparison-plots.pdf'
writer = PdfWriter()
with open(result_document, 'wb') as fh:
    writer.write(fh)

# load examples into memory
examples = []
for filepath in outpath.glob('*/'):
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

