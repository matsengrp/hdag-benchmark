"""Load tree with simulated sequences and collapse"""
import ete3

def collapse_tree(etetree, fasta):
    """Expects an ete tree with 'mutations' node attributes like that output by phastSim."""
    def node_delete(node):
        for child in node.children:
            child.mutations.extend(node.mutations)
        node.delete(prevent_nondicotomic=False)

    tree = etetree.copy()
    # collapse zero length internal edges
    for node in tree.get_descendants():
        if not node.is_leaf() and len(node.mutations) == 0:
            node_delete(node)

    # remove duplicate leaves which are grouped together
    visited = set()
    to_delete = []
    for leaf in tree.get_leaves():
        parent = leaf.up
        if id(parent) not in visited:
            visited.add(id(parent))
            identical_children = [node for node in parent.children
                                  if node.is_leaf() and len(node.mutations) == 0]
            to_delete.extend(identical_children[1:])
    for n in to_delete:
        node_delete(n)

    # remove unifurcation
    for node in tree.get_descendants():
        if len(node.children) == 1:
            node_delete(node)

    # check leaf sequences are unique
    leaves = tree.get_leaves()
    leaf_seqs = {fasta[node.name] for node in leaves}
    n = len(leaves) - len(leaf_seqs)
    if len(leaves) != len(leaf_seqs):
        raise RuntimeError("Non-unique leaf sequences in modified tree")
    return tree

def deduplicate_tree(etetree, fasta):
    """Expects an ete tree with 'mutations' node attributes like that output by phastSim.
        Same functionality as collapse_tree but keeps nodes with 0 mutations.
    """
    def node_delete(node):
        for child in node.children:
            child.mutations.extend(node.mutations)
        node.delete(prevent_nondicotomic=False)

    collapse_tree = etetree.copy()

    # collapse zero length internal edges
    for node in collapse_tree.get_descendants():
        if not node.is_leaf() and len(node.mutations) == 0:
            node_delete(node)

    # remove duplicate leaves which are grouped together
    visited = set()
    to_delete = []
    to_delete_names = []
    for leaf in collapse_tree.get_leaves():
        parent = leaf.up
        if id(parent) not in visited:
            visited.add(id(parent))
            identical_children = [node for node in parent.children
                                  if node.is_leaf() and len(node.mutations) == 0]
            to_delete.extend(identical_children[1:])
            identical_children_names = [node.name for node in parent.children
                                  if node.is_leaf() and len(node.mutations) == 0]
            to_delete_names.extend(identical_children_names[1:])
    
    for n in to_delete:
        node_delete(n)

    # remove unifurcation
    for node in collapse_tree.get_descendants():
        if len(node.children) == 1:
            node_delete(node)

    # check leaf sequences are unique
    leaves = collapse_tree.get_leaves()
    leaf_seqs = {fasta[node.name] for node in leaves}
    n = len(leaves) - len(leaf_seqs)

    if len(leaves) != len(leaf_seqs):
        raise RuntimeError("Non-unique leaf sequences in modified tree")
    
    tree = etetree.copy()
    leaves = []
    for i in range(5):
        to_delete = []
        for n in tree.get_leaves():
            if n.name in to_delete_names:
                to_delete.append(n)
        for n in to_delete:
            node_delete(n)
        
        to_delete = []
        for n in tree.get_descendants():
            if len(n.children) == 1 or (n.is_leaf() and n.name == ''):
                to_delete.append(n)
        for n in to_delete:
            node_delete(n)
        
        leaves = tree.get_leaves()
        leaf_seqs = {fasta[node.name] for node in leaves if len(node.name) > 0}
        n = len(leaves) - len(leaf_seqs)

        # print(tree)
        # print(len(leaves), len(leaf_seqs))

    if len(leaves) != len(leaf_seqs):
        raise RuntimeError("Non-unique leaf sequences in modified tree")
    
    # print(collapse_tree)
    # print("End of collapsed tree")
    # print(tree)

    return tree

def load_phastsim_newick(newickstring):
    """phastSim puts commas in node attributes, which ete can't parse"""
    tree = ete3.Tree()
    current_node = tree

    def internalnodecreate(sgen, part, current_node):
        n = current_node.add_child()
        return n

    def attributecreate(sgen, part, current_node):
        char = ''
        current_node.name = part
        while char != ']':
            _, char = next(sgen)
            part += char
        strippedlist = part.split('{')[-1].split('}')[0]
        if len(strippedlist) == 0:
            mutset = list()  # TODO: This shouldn't be a set?
        else:
            mutset = list(strippedlist.split(','))
        current_node.add_feature("mutations", mutset)
        return current_node

    def siblingcreate(sgen, part, current_node):
        pnode = current_node.up
        n = pnode.add_child()
        return n

    def moveup(sgen, part, current_node):
        return current_node.up

    def addbranchlength(sgen, part, current_node):
        idx, char = next(sgen)
        while newickstring[idx + 1] not in '([,):;':
            idx, nchar = next(sgen)
            char += nchar
        current_node.dist = float(char)
        return current_node

    tokens = {'(': internalnodecreate,
              '[': attributecreate,
              ',': siblingcreate,
              ')': moveup,
              ':': addbranchlength,
              ';': lambda _, __, cn: cn,}

    sgen = enumerate(newickstring)
    part = ''

    while True:
        try:
            current_idx, char = next(sgen)
            if char in tokens:
                current_node = tokens[char](sgen, part, current_node)
                part = ''
            else:
                part += char
        except StopIteration:
            return tree

def test_parse_newick():
    teststring = "((a,(b,c),d),(e,f));"
    teststring1 = "((a[&mutations={G28646A,C1960T}]:2.0,(b[&mutations={G28646A,C1960T}]:2.0,c[&mutations={G28646A,C1960T}]:2.0),d[&mutations={G28646A,C1960T}]:2.0)s1[&mutations={G28646A,C1960T}]:2.0,s1[&mutations={G28646A,C1960T}]:2.0(e[&mutations={G28646A,C1960T}]:2.0,f[&mutations={G28646A,C1960T}]:2.0)s1[&mutations={G28646A,C1960T}]:2.0)s1[&mutations={G28646A,C1960T}]:2.0;"
    t = ete3.Tree(teststring, format=9)
    t2 = load_phastsim_newick(teststring1)
    assert(t.write(format=9) == t2.write(format=9))

    teststring = "(s1[&mutations={G28646A,C1960T}]:2.0,((s2[&mutations={}]:1.0,s3[&mutations={G5294A,C21762T}]:1.0,(s4[&mutations={C5997T,A5570T}]:1.0,s5[&mutations={}]:0.0,s6[&mutations={}]:0.0)[&mutations={}]:1.0,s8[&mutations={G19406T,G27375T}]:1.0,s9[&mutations={C1960T}]:1.0,s10[&mutations={C3165T}]:1.0,(s11[&mutations={G20005T,G29868T}]:1.0,s12[&mutations={C5506T}]:4.0,s13[&mutations={}]:0.0,s14[&mutations={}]:0.0,s15[&mutations={}]:0.0,s16[&mutations={}]:0.0)[&mutations={}]:1.0,(s18[&mutations={}]:0.0,s19[&mutations={}]:0.0)[&mutations={C5477T,G20057T}]:1.0,(s21[&mutations={}]:0.0,s22[&mutations={}]:0.0)[&mutations={G25543A,C8917A}]:1.0,s24[&mutations={T26551C}]:1.0,s25[&mutations={G23417T}]:1.0,s26[&mutations={}]:1.0,s27[&mutations={C8172T}]:1.0,s28[&mutations={}]:0.0,s29[&mutations={}]:0.0,s30[&mutations={}]:0.0,s31[&mutations={}]:0.0,s32[&mutations={}]:0.0,s33[&mutations={}]:0.0,s34[&mutations={}]:0.0,s35[&mutations={}]:0.0,s36[&mutations={}]:0.0,s37[&mutations={}]:0.0,s38[&mutations={}]:0.0,s39[&mutations={}]:0.0,s40[&mutations={}]:0.0,s41[&mutations={}]:0.0,s42[&mutations={}]:0.0,s43[&mutations={}]:0.0,s44[&mutations={}]:0.0,s45[&mutations={}]:0.0,s46[&mutations={}]:0.0,s47[&mutations={}]:0.0,s48[&mutations={}]:0.0,s49[&mutations={}]:0.0,s50[&mutations={}]:0.0,s51[&mutations={}]:0.0,s52[&mutations={}]:0.0,s53[&mutations={}]:0.0,s54[&mutations={}]:0.0,s55[&mutations={}]:0.0,s56[&mutations={}]:0.0,s57[&mutations={}]:0.0,s58[&mutations={}]:0.0,s59[&mutations={}]:0.0,s60[&mutations={}]:0.0,s61[&mutations={}]:0.0,s62[&mutations={}]:0.0,s63[&mutations={}]:0.0,s64[&mutations={}]:0.0,s65[&mutations={}]:0.0,s66[&mutations={}]:0.0,s67[&mutations={}]:0.0,s68[&mutations={}]:0.0,s69[&mutations={}]:0.0,s70[&mutations={}]:0.0,s71[&mutations={}]:0.0,s72[&mutations={}]:0.0,s73[&mutations={}]:0.0,s74[&mutations={}]:0.0)[&mutations={}]:1.0,s76[&mutations={C3218T,G7017T}]:4.0,s77[&mutations={G27375T,G13570T}]:1.0,(s78[&mutations={}]:0.0,s79[&mutations={}]:1.0)[&mutations={C29502T}]:1.0,s81[&mutations={A10247G}]:2.0,s82[&mutations={C384T,G13403T,C29750T}]:1.0,s83[&mutations={C25240T}]:1.0,s84[&mutations={}]:1.0,(s85[&mutations={G6306T,G319T}]:5.0,s86[&mutations={}]:4.0)[&mutations={G23954T,C26845T}]:2.0,(s88[&mutations={}]:0.0,s89[&mutations={}]:0.0)[&mutations={C9447T,C17185T}]:2.0,s91[&mutations={G6884A,C321A}]:1.0,s92[&mutations={C29211T}]:1.0,s93[&mutations={A10126G,G7283T}]:1.0,s94[&mutations={G19056T}]:1.0,s95[&mutations={}]:0.0,s96[&mutations={}]:0.0,s97[&mutations={}]:0.0,s98[&mutations={}]:0.0,s99[&mutations={}]:0.0,s100[&mutations={}]:0.0,s101[&mutations={}]:0.0,s102[&mutations={}]:0.0,s103[&mutations={}]:0.0,s104[&mutations={}]:0.0,s105[&mutations={}]:0.0,s106[&mutations={}]:0.0,s107[&mutations={}]:0.0,s108[&mutations={}]:0.0,s109[&mutations={}]:0.0,s110[&mutations={}]:0.0,s111[&mutations={}]:0.0,s112[&mutations={}]:0.0,s113[&mutations={}]:0.0,s114[&mutations={}]:0.0,s115[&mutations={}]:0.0,s116[&mutations={}]:0.0,s117[&mutations={}]:0.0,s118[&mutations={}]:0.0,s119[&mutations={}]:0.0,s120[&mutations={}]:0.0,s121[&mutations={}]:0.0,s122[&mutations={}]:0.0,s123[&mutations={}]:0.0,s124[&mutations={}]:0.0,s125[&mutations={}]:0.0,s126[&mutations={}]:0.0,s127[&mutations={}]:0.0,s128[&mutations={}]:0.0)[&mutations={G21424A,G10051A,C27845T}]:1.0,s130[&mutations={}]:1.0)[&mutations={G862T,G16390T,G10739T,C26343T,T20746A,G669T}]:7.0;\n"

    teststring1 = ''
    g = iter(teststring)
    while True:
        try:
            char = next(g)
            if char == '[':
                while char != ']':
                    char = next(g)
            else:
                teststring1 += char
        except StopIteration:
            break
    t = ete3.Tree(teststring1)
    t2 = load_phastsim_newick(teststring)
    assert(t.write(format=9) == t2.write(format=9))
    return t2

    
