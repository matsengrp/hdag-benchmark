import click
import historydag as hdag
import ete3 as ete
import random
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.metrics import precision_recall_curve
import numpy as np


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Command line scripts to facilitate node support validation
    """
    pass


###################################################################################################
### Simulation Info ###############################################################################
###################################################################################################

# TODO: Use this in all simulations scripts
@click.command("get_pars_score")
@click.option('--sim_dir', '-s', help='the folder containing TOI and fasta file.')
def get_pars_score(sim_dir):
    """
    Load true tree as an ete3 tree called `tree`, with leaf sequences
    stored in nodes' attribute `sequence`
    """
    from historydag.parsimony import parsimony_score, sankoff_upward
    tree_path = sim_dir + "/collapsed_simulated_tree.nwk"
    tree = ete.Tree(tree_path) # Doesn't have internal names

    fasta_path = sim_dir + "/ctree_with_refseq.fasta"   # ancestral seq in second line of this file
    with open(fasta_path, "r") as f:
        assert ">ancestral\n" == f.readline()
        ancestral_seq = f.readline().strip()
    
    # build sequences from mutations
    for node in tree.traverse("preorder"):
        if node.is_root():
            seq = ancestral_seq
        else:
            seq = node.up.sequence
            if len(node.mutations) >= 1:
                for mut in node.mutations.split("|"):
                    mut = mut.strip()
                    curr = mut[0]
                    i = int(mut[1:-1])-1    # convert indices to 0-based
                    new = mut[-1]
                    assert seq[i] == curr
                    seq = seq[:i] + new + seq[i + 1:]
        
        node.add_feature("sequence", seq)

    # just counts mutations between simulated internal node sequences
    tree_score = parsimony_score(tree)

    # computes the best possible parsimony score of any labeling on tree's topology
    max_score = sankoff_upward(tree, len(tree.sequence))

    outfile = sim_dir + "/toi_info/tree_stats.txt"
    with open(outfile, "w") as f:
        f.write(f"TOI parsimony score\t{tree_score}\n")
        f.write(f"Max parsimony score (topology)\t{max_score}\n",)
        f.write(f"Max parsimony score (on data)\t...Use dnapars script\n")




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
        
        etenode2cu[node.name] = frozenset(cu)

    return set([v for k, v in etenode2cu.items()])

def hdag_output(node_set, pb_file, taxId2seq):
    """
    Returns a list of tuples (clade_id, estimated_support, in_tree). The list is primarily
    sorted by estimated support, and portions that have the same support are randomly shuffled
    """
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(pb_file)
    # TODO: Comment out this line... maybe
    dag.make_complete()
    dag.trim_optimal_weight() # Trim to the MP trees
    counts = dag.count_nodes(collapse=True) # TODO: double check that this is counting the right thing
    total_trees = dag.count_trees()
    

    print("size of counts is:", len(counts))   # TODO: Check how many nodes are in the counts
    print(f"There are {len(node_set)} nodes in TOI")
    dag.summary()

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
        # TODO: Fix key error labels seq is not among taxID
        id_node = frozenset([seq2taxId[label.compact_genome.to_sequence()] for label in node])
        est_sup = counts[node] / total_trees
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

import subprocess


# TODO: Implement a way to detect if you are slowing down, and then add the sample from any tree option    
#
# $larch_usher_exec -i $seedtree -r $refseqfile -c 400 -o $optdag1 --move-coeff-pscore 1 --move-coeff-nodes 1 -l $log1
@click.command("larch_usher")
@click.option('--executable', '-e', default='/home/whowards/larch/larch/build/larch-usher', help='path to pre-made larch-usher executable')
@click.option('--input', '-i', help='input tree or hdag. if tree, need refseqfile.')
@click.option('--refseqfile', '-r', default=None, help='number of .')
@click.option('--count', '-c', help='number of iterations.')
@click.option('--out_dir', '-o', help='the directory for where to store resulting dag protobufs.')
@click.option('--schedule', '-s', default="annealed")
@click.option('--log_dir', '-l')
@click.option('--pars_score', '-p', default=1)
@click.option('--node_score', '-n', default=1)
def larch_usher(executable, input, refseqfile, count, out_dir, schedule, log_dir, pars_score, node_score):
    """Python CLI for driving larch-usher"""

    os.chdir(f"{out_dir}")
    # subprocess.run(["cd", out_dir]) # One process can't change anothers working dir

    if int(count) <= 2:
        return

    # Cast a wide net by prioritizing new nodes only
    print("Running initial iterations of larch-usher...")
    subprocess.run(["mkdir", "-p", f"{log_dir}_1"])
    args = [executable,
            "-i", input,
            "-c", f"{round(int(count)/2)}",
            "-o", f"{out_dir}/opt_dag_1.pb",
            "-l", f"{log_dir}_1",
            "--move-coeff-nodes", str(1),
            "--move-coeff-pscore", str(0),
            "--sample-best-tree"            # NOTE: Might need to change this with different version of larch-usher
            ]
    if refseqfile is not None:
        args.extend(["-r", refseqfile])
    subprocess.run(args=args)

    # Start considering parsimonious moves
    subprocess.run(["mkdir", "-p", f"{log_dir}_2"])
    args = [executable,
            "-i", f"{out_dir}/opt_dag_1.pb",
            "-c", f"{round(int(count)/6)}",
            "-o", f"{out_dir}/opt_dag_2.pb",
            "-l", f"{log_dir}_2",
            "--move-coeff-nodes", str(2),
            "--move-coeff-pscore", str(1),
            "--sample-best-tree"
            ]
    subprocess.run(args=args)

    # Emphasize parsimony over new nodes
    subprocess.run(["mkdir", "-p", f"{log_dir}_3"])
    args = [executable,
            "-i", f"{out_dir}/opt_dag_2.pb",
            "-c", f"{round(int(count)/6)}",
            "-o", f"{out_dir}/opt_dag_3.pb",
            "-l", f"{log_dir}_3",
            "--move-coeff-nodes", str(1),
            "--move-coeff-pscore", str(3),
            # "--sample-any-tree"
            ]
    subprocess.run(args=args)

    print("Completing DAG...")
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(f"{out_dir}/opt_dag_3.pb")
    dag.make_complete()
    dag.to_protobuf_file(f"{out_dir}/complete_opt_dag.pb")

    subprocess.run(["mkdir", "-p" f"{log_dir}_complete"])
    args = [executable,
            "-i", f"{out_dir}/complete_opt_dag.pb",
            "-c", f"{round(int(count)/6)}",
            "-o", f"{out_dir}/final_opt_dag.pb",
            "-l", f"{log_dir}_complete",
            "--move-coeff-nodes", str(1),
            "--move-coeff-pscore", str(3),
            "--sample-best-tree"
            ]
    subprocess.run(args=args)


    # TODO: Add --sample-any again?




def beast_output(node_sets, pb_file, seq2taxId):
    """... """
    return None


###################################################################################################
#### Aggregation ##################################################################################
###################################################################################################

# TODO: Rename this
@click.command("agg")
@click.option('--input', '-i', help='file path to input list as a pickle.')
@click.option('--out_dir', '-o', help='output directory to store figures/tables in.')
@click.option('--clade_name', '-c')
def agg(input, out_dir, clade_name):
    """Given the pickled file, aggregates results for support values"""

    with open(input, "rb") as f:
        results = pickle.load(f)

    window_size = int(len(results) * 0.20)
    out_path = out_dir + f"/fig_w={window_size}.png"
    x, y = sliding_window_plot(results, window_size=window_size)

    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[0.7, 0.3])
    f.set_size_inches(7, 9)
    f.suptitle(f"Clade {clade_name} Coverage Analysis")

    ax1.set_ylabel(f"Empirical Probability (window_size={window_size}/{len(results)})")
    ax1.plot(x, y)
    ax1.plot([0, 1], [0, 1])
    
    ax2.hist(x)
    ax2.set_xlabel("Estimated Support")
    # ax2.set_yscale("log")
    f.savefig(out_path)
    f.clf()

    # TODO: Plot PR curves in different method
    leaf_idxs = []
    for i, (clade, _, _) in enumerate(results):
        if len(clade) <= 1:
            leaf_idxs.append(i)
    

    f.set_size_inches(6.4, 4.8)
    y_true = [y for i, (_, _, y) in enumerate(results)]
    y_scores = [sup for i, (_, sup, _) in enumerate(results)]
    ratio = sum(y_true) / len(y_true)

    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
    # print(precision)
    # print(recall)
    # print(thresholds)
    plt.ylabel("Precision (TP / TP + FP)")      # intuitively the ability of the classifier not to label as positive a sample that is negative.
    plt.xlabel("Recall (TP / TP + FN)")         # intuitively the ability of the classifier to find all the positive samples
    plt.plot(recall, precision, label="Support Classifier")
    plt.plot([0, 1], [ratio, ratio], label="Random")             # Random classifier performance (i.e., proportion of positives)
    plt.legend()
    plt.title(f"PR Curve")
    plt.xlim(0, 1.02)
    plt.ylim(0, 1.02)
    plt.savefig(out_dir + f"/prec_rec.png")
    plt.clf()

    # TODO: Plot threshold curve
    leaf_idxs = []
    for i, (clade, _, _) in enumerate(results):
        if len(clade) <= 1:
            leaf_idxs.append(i)
    

    f.set_size_inches(6.4, 4.8)
    y_true = [y for i, (_, _, y) in enumerate(results)]
    y_scores = [sup for i, (_, sup, _) in enumerate(results)]
    ratio = sum(y_true) / len(y_true)


    thresholds = np.linspace(0, 1, 21)
    acc_sup = []
    for thresh in thresholds:
        correct = 0
        for y, sup in zip(y_true, y_scores):
            if (sup >= thresh and y == 1) or (sup < thresh and y == 0):
                correct+= 1
        
        acc_sup.append(correct/ len(y_true))

    plt.ylabel("Accuracy")
    plt.xlabel("Threshold")
    plt.plot(thresholds, acc_sup, label="Support Classifier")
    plt.plot([0, 1], [ratio, ratio], label="Random")             # Random classifier performance (i.e., proportion of positives)
    plt.legend()
    plt.title(f"Accuracy vs Threshold")
    plt.xlim(0, 1.02)
    plt.ylim(ratio-0.02, 1.02)
    plt.savefig(out_dir + f"/thresh_eval.png")







@click.command("clade_results")
@click.option('--clade_dir', '-c', help='path to clade directory.')
@click.option('--num_sim', '-n', default=1, help='number of simulations to average.')
def clade_results(clade_dir, num_sim, method="historydag"):
    """Given the clade directory, performs coverage analysis across all simulations"""

    clade_name = clade_dir.split("/")[-1]

    result_dict = {}
    for trial in range(1, num_sim+1):
        # Assumes that `path/to/clade/trial/results/results.pkl stores`` list of nodes and their supports
        result_path = clade_dir + f"/{trial}/results/historydag/results.pkl"
        with open(result_path, "rb") as f:
            results = pickle.load(f)
        result_dict[trial] = results

    window_size = int(len(results) * 0.2)
    out_path = clade_dir + f"/figures/multi_line_w={window_size}.png"

    for trial in range(1, num_sim+1):
        x, y = sliding_window_plot(result_dict[trial], window_size=window_size)
        plt.plot(x, y)
    plt.plot([0, 1], [0, 1])
    plt.xlabel("Estimated Support")
    plt.ylabel(f"Empirical Probability (window_size={window_size}/{len(results)})")
    plt.title(f"Aggregated {clade_name} Coverage Analysis")
    plt.savefig(out_path)
    plt.clf()


    
    out_path = clade_dir + f"/figures/single_line_w={window_size}.png"

    results_full = []
    for k, v in result_dict.items():
        results_full.extend(v)
    
    results_full.sort(key=lambda el: el[1])
    x, y, pos_devs, neg_devs = sliding_window_plot(results_full, std_dev=True, window_size=window_size)
    plt.plot(x, y)
    plt.plot(x, pos_devs, linestyle='dashed', color='green')
    plt.plot(x, neg_devs, linestyle='dashed', color='green')
    plt.plot([0, 1], [0, 1])
    plt.xlabel("Estimated Support")
    plt.ylabel(f"Empirical Probability (window_size={window_size}/{len(results)})")
    plt.title(f"Aggregated {clade_name} Coverage Analysis")
    plt.savefig(out_path)


def sliding_window_plot(results, std_dev=False, window_size=200):
    """Given list of results tuples returns xy coords of sliding window plot."""

    x, y = [], []
    devs = []
    side_len = int(window_size/2)
    for i, (_, est_sup, _) in enumerate(results):
        x.append(est_sup)
        window = [int(el[2]) for el in results[max(0, i-side_len):min(len(results), i+side_len)]]
        y.append(sum(window) / len(window))

        if std_dev:
            avg = y[i]
            sq_diff = [(el - avg)**2 for el in window]
            devs.append((sum(sq_diff) / len(sq_diff))) #** 0.5) # TODO: Should we be using variance or std dev
    
    if std_dev:
        pos_devs = [y_val + dev for y_val, dev in zip(y, devs)]
        neg_devs = [y_val - dev for y_val, dev in zip(y, devs)]
        return x, y, pos_devs, neg_devs
    else:
        return x, y



cli.add_command(get_pars_score)
cli.add_command(clade_results)
cli.add_command(save_supports)
cli.add_command(larch_usher)
cli.add_command(agg)

if __name__ == '__main__':
    cli()