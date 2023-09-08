from collections import Counter
import historydag as hdag

def get_muts(node):
    if node.is_ua_node():
        yield from ()
    else:
        pnode = next(iter(node.parents))
        if pnode.is_ua_node():
            yield from ()
        else:
            if "sequence" in pnode.label._fields:
                for idx, (pbase, cbase) in enumerate(zip(pnode.label.sequence, node.label.sequence)):
                    if pbase != cbase:
                        yield pbase + str(idx + 1) + cbase
            else:
                for pbase, cbase, site in hdag.compact_genome.cg_diff(pnode.label.compact_genome, node.label.compact_genome):
                    yield pbase + str(site) + cbase

def predict_hotspots(dag):
    site_counter = Counter()
    mut_counter = Counter()
    for node in dag.preorder(skip_ua_node=True):
        muts = list(get_muts(node))
        site_counter.update(int(mut[1:-1]) for mut in muts)
        mut_counter.update(muts)
    return site_counter, mut_counter
