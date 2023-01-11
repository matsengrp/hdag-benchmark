import click
import historydag as hdag
import ete3 as ete
import random
import pickle
import matplotlib.pyplot as plt
import seaborn as sns


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Command line scripts to facilitate node support validation
    """
    pass

###################################################################################################
#### Inference ####################################################################################
###################################################################################################

@click.command("save_supports")
@click.option('--method', '-m', default="hdag", help='method of inference.')
@click.option('--tree_path', '-t', help='the newick file for the true tree.')
@click.option('--input_path', '-i', help='the file containing inference results.')
@click.option('--output_path', '-o', help='the file to save output to.')
def save_supports(method, tree_path, input_path, output_path):
    """
    A method for computing, formatting, and saving node supports for various methods of inference 
    Input:  method of inference (e.g., hdag, beast, etc.),
            input file path (e.g., to history dag, or beast output),
            output filepath 
    Output: For each file in the input directory, a text file containing list of support values in the
            format (clade, estimated_support, in_tree)
    """

    # Map of taxon id (e.g., s1, s4, etc) to full sequence
    fasta_path = tree_path + ".fasta"
    taxId2seq = hdag.parsimony.load_fasta(fasta_path)

    # Compute the set of nodes that are in the true tree
    node_set = get_true_nodes(tree_path)

    # Computes list of nodes 
    if method == "hdag":
        support_list = hdag_output(node_set, input_path, taxId2seq)
    elif method == "beast":
        support_list = beast_output(node_set, input_path, seq2taxId)
    else:
        assert False # TODO: Throw error

    with open(output_path, "wb") as f:
        pickle.dump(support_list, f)
    

def get_true_nodes(tree_path):
    """ Given a newick file, returns the set of all clades (as frozen sets of taxon ids)
    """
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
        
        print(cu)
        
        etenode2cu[node.name] = frozenset(cu)

    return set([v for k, v in etenode2cu.items()])

def hdag_output(node_set, pb_file, taxId2seq):
    """
    Returns a list of tuples (clade_id, estimated_support, in_tree). The list is primarily
    sorted by estimated support, and portions that have the same support are randomly shuffled
    """
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(pb_file)
    dag.trim_optimal_weight() # Trim to the MP trees
    counts = dag.count_nodes(collapse=True) # TODO: double check that this is counting the right thing
    total_trees = dag.count_trees()

    # dag.summary()
    print(f"Dag contains {total_trees} MP trees")
    print("size of counts is:", len(counts))   # TODO: Check how many nodes are in the counts
    print(f"There are {len(node_set)} nodes in TOI")

    seq2taxId = {v: k for k, v in taxId2seq.items()}

    # Construct results dict that maps nodes (frozen sets of taxon ids) to tuples of estimated
    # support and whether that node is in the true tree or not
    node2stats = {}

    # Get the support for all dag nodes
    counter = 0
    for node in counts:
        if len(node) == 0:  # UA node
            continue

        # Convert cg label to taxon id
        id_node = frozenset([seq2taxId[label.compact_genome.to_sequence()] for label in node])
        est_sup = counts[node] / total_trees
        node2stats[id_node] = (est_sup, id_node in node_set)
    

    # Get the support for all nodes in true tree
    for id_node in node_set:
        if id_node not in node2stats.keys():
            node2stats[id_node] = (0, True)

    print("Coinsidering", len(node2stats), "nodes")
    stats_list =[(id_node, stats[0], stats[1]) for id_node, stats in node2stats.items()]
    random.shuffle(stats_list)
    stats_list.sort(key=lambda el: el[1])

    return stats_list


def beast_output(node_sets, pb_file, seq2taxId):
    """... """
    return None


###################################################################################################
#### Aggregation ##################################################################################
###################################################################################################

@click.command("agg")
@click.option('--input', '-i', help='file path to input list as a pickle.')
@click.option('--out_dir', '-o', help='output directory to store figures/tables in.')
def agg(input, out_dir):
    """Given the pickled file, aggregates results for support values"""
   
    with open(input, "rb") as f:
        results = pickle.load(f)

    counter = 0
    for node, sup, in_tree in results:
        if len(node) > 1:
            print(sup, "\t" ,in_tree, "\t", node)
            counter += 1
    print(counter, "nodes under consideration")

    out_path = out_dir + "/fig.png"
    x, y = sliding_window(results, window_size=200)
    plt.plot(x, y)
    plt.plot([0, 1], [0, 1])
    plt.savefig(out_path)


def sliding_window(results, window_size=20):
    """Given list of results tuples returns xy coords of sliding window plot."""

    x, y = [], []
    side_len = int(window_size/2)
    for i, (_, est_sup, _) in enumerate(results[side_len+1 : -side_len-1]):
        idx = side_len+1+i  # index of starting position in list
        x.append(est_sup)
        window = [int(el[2]) for el in results[idx-side_len:idx+side_len]]
        print(idx, window)
        y.append(sum(window) / len(window))
    
    return x, y

    



cli.add_command(save_supports)
cli.add_command(agg)

if __name__ == '__main__':
    cli()