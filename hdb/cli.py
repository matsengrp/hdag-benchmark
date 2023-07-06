"""Command line interface."""

import sys
import click
import hdb.summary as summary
import hdb.collapse_tree as ct
import hdb.aggregate_dnapars_trees as agg
import ete3
import pickle
import random
import math

# NOTE: These should probably be moved to a different folder..?
import subprocess
from collections import Counter
import ete3 as ete
from historydag.parsimony import parsimony_score, sankoff_upward
import json

# Entry point
def safe_cli():
    """Top-level CLI for subcommands."""
    try:
        cli()
    except Exception as exception:
        print("Exception raised when running the command:\n")
        print(" ".join(sys.argv) + "\n")
        raise exception


@click.group(
    context_settings={"help_option_names": ["-h", "--help"]},
    invoke_without_command=True,
)
def cli():
    """Top-level CLI for subcommands."""
    pass  # pylint: disable=unnecessary-pass


@cli.command()
@click.argument("outfiles", nargs=-1)
@click.option("-r", "--root", help="name for outgroup sequence")
@click.option("-a", "--abundance_file", help="filepath to abundance mapping")
@click.option("-o", "--output_path", help="filepath to write pickled hDAG")
def aggregate_dnapars_trees(outfiles, root, abundance_file, output_path):
    """Process dnapars outfiles into a history DAG."""
    dag = agg.aggregate_dnapars_trees(outfiles, root, abundance_file)
    with open(output_path, 'wb') as fh:
        pickle.dump(dag, file=fh)

@cli.command()
@click.argument("input_path")
@click.argument("out_path_base")
def numerify_taxon_names(input_path, out_path_base):
    """Convert taxon names to s<number>, outputting a "mapping file" and a renamed newick."""

    mapping_path = out_path_base + ".mapping"
    numerified_path = out_path_base + ".n.nwk"

    tree = ete3.Tree(input_path, format=1)
    mapping = {}
    curr_num = 1
    for node in tree.traverse("postorder"):
        mapping[node.name] = f"s{curr_num}"
        node.name = f"s{curr_num}"
        curr_num += 1

    tree.write(format=1, outfile=numerified_path)
    with open(mapping_path, "w") as f:
        for seq_id, tax_id in mapping.items():
            print(f"{seq_id} {tax_id}", file=f)

@cli.command()
@click.argument("indags", nargs=-1)
@click.option("-o", "--output_path", help="filepath to write pickled hDAG")
def merge_dags(indags, output_path):
    def load_dag(indag):
        with open(indag, 'rb') as fh:
            dag = pickle.load(fh)
        return dag
    
    dag = agg.merge_dags(load_dag(indag) for indag in indags)
    with open(output_path, 'wb') as fh:
        pickle.dump(dag, file=fh)

@cli.command()
@click.argument("fasta_paths", nargs=-1)
@click.option("-o", "--output-path", default="summary.csv", help="CSV output path.")
def alnsummarize(fasta_paths, output_path):
    """Summarize a collection of alignments via AMAS."""
    df = summary.summary_df_of_fasta_paths(fasta_paths)
    df.to_csv(output_path)

@cli.command()
@click.argument("input_newick")
@click.argument("input_fasta")
@click.argument("output_newick")
def collapse_tree(input_newick, input_fasta, output_newick):
    with open(input_newick, 'r') as fh:
        intree = ct.load_phastsim_newick(fh.read())
    infasta = load_fasta(input_fasta)
    outtree = ct.collapse_tree(intree, infasta)
    # # TODO:
    # ct.get_tree_stats(sim_dir)

    outtree.write(features=["mutations"], format_root_node=True, outfile=output_newick)
    variant_sites = set()
    for node in outtree.traverse():
        for mut in node.mutations:
            variant_sites.add(int(mut[1:-1]))
    variant_sites = list(sorted(variant_sites))
    with open(output_newick + '.variant_sites.txt', 'w') as fh:
        fh.write(' '.join(str(num) for num in variant_sites))
    
    def excise_variants(seq):
        return ''.join(seq[idx - 1] for idx in variant_sites)

    with open(output_newick + '.fasta', 'w') as fh, open(output_newick + '.variant_sites.fasta', 'w') as fhvariants:
        for seqname in sorted((n.name for n in outtree.get_leaves()),
                              key=lambda name: int(name[1:])):
            print('>' + seqname, file=fh)
            print('>' + seqname, file=fhvariants)
            print(infasta[seqname], file=fh)
            print(excise_variants(infasta[seqname]), file=fhvariants)

# TODO: Maybe this should go in ct
@cli.command()
@click.option('--sim_dir', '-s', help='the folder containing TOI and fasta file.')
def get_tree_stats(sim_dir):
    """
    Computes the parsimony score of the simulated tree, the maximum possible parsimony given the
    topology, and the maximum parsimony on the leaves. Stores results as a json at sim_dir/tree_stats.json

    Also computes statistics for the subsetted USHER tree. This can be useful for ensuring that
    the simulations roughly match the real data. This script assumes the directory `sim_dir/..`
    contains the USHER-subsetted newick tree tree.n.nwk.
    """

    var_sites_prefix = sim_dir + "/collapsed_simulated_tree.nwk.variant_sites"
    with open(var_sites_prefix + ".txt", "r") as f:
        idxs = f.readline().split()

    with open(sim_dir + "/ctree_with_refseq.fasta", "r") as f:
        name = f.readline().strip()[1:] # NOTE: Should be `ancestral`
        seq = f.readline()              # Entire sequence
        variants = ""                   # The character values of the nucs that vary
        for i in idxs:
            variants += seq[int(i)-1]

    with open(var_sites_prefix + ".fasta", "r") as f:
        var_sites = f.readlines()

    with open(var_sites_prefix + "_with_refseq.fasta", "w") as f:
        f.write(f">{name}\n")
        f.write(f"{variants}\n")
        for line in var_sites:
            f.write(f"{line}")

    # subprocess.run([
    #     "dnapars_parsimony_score.sh",
    #     var_sites_prefix + "_with_refseq.fasta",    # infasta
    #     name,                                       # root_name
    #     sim_dir                                     # out_dir
    # ])
    # with open(sim_dir + "/dnapars_output.txt", "r") as f:
    #     line = f.readline()
    #     temp = line.strip().split(" ")
    #     best_possible = float(temp[-1])
    best_possible = -1
    
    tree_path = sim_dir + "/collapsed_simulated_tree.nwk"
    tree = ete.Tree(tree_path) # Doesn't have internal names

    multifurc_counts = Counter()
    for node in tree.traverse():
        if not node.is_leaf():
            multifurc_counts[len(node.children)] += 1

    fasta_path = sim_dir + "/ctree_with_refseq.fasta"   # ancestral seq in second line of this file
    with open(fasta_path, "r") as f:
        assert ">ancestral\n" == f.readline()
        ancestral_seq = f.readline().strip()
    
    # build sequences from mutations
    num_nodes = 0
    num_leaves = 0
    for node in tree.traverse("preorder"):
        num_nodes += 1
        if node.is_leaf():
            num_leaves += 1
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

    # computes the best possible parsimony score of any labeling on tree's
    # topology, with root sequence constrained.
    # first replace all internal, non-root sequences with N's, then do
    # sankoff_upward with use_internal_node_sequences True.
    sankoff_tree = tree.copy()
    for node in sankoff_tree.traverse():
        if node.is_root() or node.is_leaf():
            continue
        else:
            node.sequence = "N" * len(tree.sequence)
    max_score = sankoff_upward(sankoff_tree, len(sankoff_tree.sequence), use_internal_node_sequences=True)

    stats_dict = {
        "num_leaves": num_leaves,
        "num_nodes": num_nodes,
        "pars_score": tree_score,
        "max_score_top": max_score,
        "max_score_data": best_possible,
        "multifurc_distribution": multifurc_counts
    }

    print(f"MP with topology: {max_score}\t MP on data: {best_possible}")

    # Add similar stats from the USHER tree
    clade_dir = sim_dir.split("/")[0]
    usher_tree_path = clade_dir + "/tree.n.nwk"
    usher_tree = ete.Tree(usher_tree_path, format=1)
    # Remove multifurcations
    to_delete = []
    for node in usher_tree.traverse():
        # TODO: Add leaf check elsewhere as needed
        if not node.is_root() and node.dist==0:
            to_delete.append(node)
    for node in to_delete:
        node.delete(prevent_nondicotomic=False)
    # Remove unifurcations
    to_delete = [node for node in usher_tree.traverse() if len(node.children) == 1 and not node.is_root()]
    for node in to_delete:
        node.delete(prevent_nondicotomic=False)

    pars = 0
    multifurc_counts = Counter()
    node_count = 0
    leaf_count = 0
    for node in usher_tree.traverse():
        node_count += 1
        if node.is_leaf():
            leaf_count += 1
        if not node.is_root():
            pars += node.dist
        if not node.is_leaf():
            multifurc_counts[len(node.children)] += 1
    
    stats_dict["usher_num_leaves"] = leaf_count
    stats_dict["usher_num_nodes"] = node_count
    stats_dict["usher_pars_score"] = pars
    stats_dict["usher_multifurc_distribution"] = multifurc_counts



    outfile = sim_dir + "/tree_stats.json"
    with open(outfile, "w") as f:
        f.write(json.dumps(stats_dict, indent=4))

    # TODO: This is wrong, right?
    # if best_possible != max_score:
    #     raise RuntimeError("Non-unique leaf sequences in modified tree")


@cli.command()
@click.option("-i", "--input-path", help="Newick tree input path.")
@click.option("-o", "--output-path", help="Output path.")
@click.option("-s", "--resolve-seed", default=1)
@click.option("-b", "--branch-len-model", default="num-muts")
@click.option('--add-ancestral/--no-add-ancestral', default=False, help='add leaf below root labeled `ancestral`')
def resolve_multifurcations(input_path, output_path, resolve_seed, branch_len_model, add_ancestral):
    """
    Given an ete tree, resolves all polytomies by creating a
    uniformly random bifurcating tree that is consistent with
    the multifurcating one.
    """
    tree = ete3.Tree(input_path, format=1)
    resolve_polytomy(tree, resolve_seed, branch_len_model)
    if add_ancestral:
        new_tree = ete3.Tree()
        new_tree.add_child(name='ancestral', dist=0)
        new_tree.add_child(child=tree, dist=0)
        tree = new_tree
    tree.write(outfile=output_path, format=1)


def resolve_polytomy(tree, seed, branch_len_model):
    """
    Given an ete tree, resolves all polytomies by creating a
    uniformly random bifurcating tree that is consistent with
    the multifurcating one.
    """

    # Controls how much multifurcation in simulation
    # TODO: Trying smaller values here (was 0.01)
    resolved_multifurc_len = 0.001

    with open("refseq.fasta", "r") as f:
        f.readline()
        genome = f.readline()
        genome_len = len(genome)
        print("Genome has length", genome_len)

    random.seed(seed)
    new_node_name = 1

    def jc_dist(num_mut):
        # genome_len = 29904
        mut_prop = num_mut / genome_len
        return -3/4 * math.log(1 - 4/3 * mut_prop)

    def _resolve(node, new_node_name):
        if len(node.children) > 2:
            node_list = list(node.children)
            avg_mut = sum([n.dist for n in node_list]) / len(node_list) # Average number of mutations from parent
            node.children = []
            while len(node_list) > 2:
                # Randomly sample a pair of nodes
                pair = random.sample(range(0, len(node_list)-1), 2)
                pair.sort() # TODO: Why are we doing this?
                
                # Create new node and make branch length
                par = ete3.Tree()

                if branch_len_model == "jc":
                    adj_dist = resolved_multifurc_len
                    par.dist = jc_dist(resolved_multifurc_len)
                elif branch_len_model == "num-muts":
                    adj_dist = resolved_multifurc_len
                    par.dist = adj_dist

                # merge under a parent node
                par.name = f"r{new_node_name}"
                new_node_name += 1
                par.add_child(node_list.pop(pair[1]))
                par.add_child(node_list.pop(pair[0]))

                # insert into list for this node to be further merged
                node_list.append(par)
            
            node.add_child(node_list[0])
            node.add_child(node_list[1])

    target = [tree]
    # TODO: Should we use get_descendants instead of postorder??
    # target.extend([n for n in tree.traverse("postorder")])
    target.extend([n for n in tree.get_descendants()])
    for n in target:
        # Compute branch lengths for edge above node if not a parent node
        if not n.is_root():
            if branch_len_model == "jc":
                branch_len = jc_dist(n.dist)
            elif branch_len_model == "num-muts":
                branch_len = n.dist

            n.dist = branch_len
        
        _resolve(n, new_node_name) # Resolve the edges under given node



def load_fasta(fastapath):
    fasta_records = []
    current_seq = ''
    with open(fastapath, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                fasta_records.append([line[1:].strip(), ''])
            else:
                fasta_records[-1][-1] += line.strip()
    return dict(fasta_records)

# @cli.command()
# @click.option("-c", "--clade-path", help="Path to clade directory.")
# @click.option("-n", "--num-trials", default=25, help="Number of trials to return.")
# def find_diversity(clade_path, num_trials):
#     trial_list = []
#     for trial in range(1, 101):
#         log_path = clade_path + f"/{trial}/results/historydag/opt_info/optimization_log_complete/logfile.csv"
#         try:
#             with open(log_path, "r") as f:
#                 num_trees = f.readlines()[-1].split("\t")[-3]
#                 trial_list.append((trial, int(num_trees)))
#         except:
#             print(f"\tSkipping {log_path}")
#             continue
    
#     trial_list.sort(key=lambda el: el[1], reverse=True)
#     trial_list = trial_list[:num_trials]
#     with open("pars_div_trials.txt", "w") as f:
#         for trial, num_trees in trial_list:
#             f.write(f"{trial}\n")
#             print(trial, "\t", num_trees)
            

@cli.command()
@click.option("-c", "--clade-path", help="Path to clade directory.")
@click.option("-n", "--num-trials", default=25, help="Number of trials to return.")
def get_trial_info(clade_path, num_trials):
    trial_list = []
    for trial in range(1, 26):
        log_path = clade_path + f"/{trial}/results/historydag/opt_info/optimization_log_complete/logfile.csv"   # NOTE: This is somewhat unreliable measure of parsimony diversity
        sim_info_path = clade_path + f"/{trial}/simulation/tree_stats.json"
        results_path = clade_path + f"/{trial}/results/historydag/results.pkl"
        support_log_path = clade_path + f"/{trial}/results/inference__historydag.log"
        dnapars_path = clade_path + f"/{trial}/simulation/dnapars/outfile"
        
        try:
            with open(results_path, "rb") as f:
                results = pickle.load(f)
                results_in_hdag = [res for res in results if res[1] > 0]
                num_nodes = len(results_in_hdag)

            with open(sim_info_path, "r") as f:
                tree_stats = json.load(f)
                toi_score = int(tree_stats['max_score_top'])

            with open(log_path, "r") as f:
                last_line = f.readlines()[-1]
                num_trees_uncollapsed = int(last_line.split("\t")[-3])
                max_dag_pars = int(last_line.split("\t")[-4])

            with open(support_log_path, "r") as f:
                line = f.read().split("=>")[-1].strip()
                assert "DAG contains " in line
                num_trees = int(line.split(" ")[2])
                counter_text = line.split("\n")[3][9:-2]
                kv_text_list = counter_text.split(", ")
                dists = {}
                for string in kv_text_list:
                    num_list = string.split(": ")
                    dist = int(num_list[0])
                    count = int(num_list[1])
                    dists[dist] = count
                closest_dist = min(list(dists.keys()))
                max_dist = max(list(dists.keys()))

            with open(dnapars_path, "r") as f:
                line = f.readlines()[4].strip()
                dnapars_trees = int(line.split(" ")[0])
            
            # List of trial number, number of trees in uncollapsed DAG, toi and max pars scores, number of nodes in final DAG
            if closest_dist >= 0:
                trial_list.append((trial, num_trees_uncollapsed, num_trees, toi_score, max_dag_pars, num_nodes, closest_dist, max_dist, dnapars_trees))
        except Exception as e:
            print(e)
            continue

    trial_list.sort(key=lambda el: el[-3], reverse=True)

    if num_trials < len(trial_list):
        trial_list = trial_list[:num_trials]

    with open(f"{clade_path}/pars_div_trials.txt", "w") as f:
        print('trial, num_trees, toi_score, max_dag_pars, num_nodes, closest_dist, max_dist, dnapars_trees')
        for trial, num_trees_uncollapsed, num_trees, toi_score, max_dag_pars, num_nodes, closest_dist, max_dist, dnapars_trees in trial_list:
            f.write(f"{trial}\n")
            print(trial, "\t", \
                    # num_trees_uncollapsed, "\t", \
                    num_trees, "\t", \
                    toi_score, "\t", \
                    max_dag_pars, "\t", \
                    num_nodes, "\t", \
                    closest_dist, "\t", \
                    max_dist, "\t", \
                    dnapars_trees
                )

# @cli.command()
# @click.option("-c", "--clade-path", help="Path to clade directory.")
# def reconstruct_fasta(clade_path):
#     """
#     Reconstructs 'original' fasta file from tree and leaf sequences. Expects USHER protobuf to
#     be stored at `clade_tree.pb.gz` and the reference sequence to be at `refseq.fasta` both in the
#     `clade_path` directory.

#     NEED to run `conda activate bte` first
#     """

#     import bte
#     tree = bte.MATree(f'{clade_path}/clade_tree.pb.gz')

#     # Build new reference sequence
#     with open(f'{clade_path}/refseq.fasta', 'r') as f:
#         f.readline()
#         refseq_orig = f.readline()
#     mutations = {}
#     for node in tree.breadth_first_expansion():
#         for mut in node.mutations:
#             base = mut[0]
#             idx = int(mut[1:-1])-1
#             if idx not in mutations:
#                 mutations[idx] = base
#     refseq = edit_string(refseq_orig, mutations)

#     # Create edit dictionary for each node and apply edits
#     leaf2seq = {}
#     leaves = tree.get_leaves()
#     for leaf in leaves:
#         mutations = []
#         mutations.extend(leaf.mutations)
#         node = leaf
#         while node.level > 1:
#             node = node.parent
#             mutations.extend(node.mutations)
#             if node.level == 1:
#                 break
        
#         mutations.reverse()
#         node_seq = refseq
#         muts = {}
#         for mut in mutations:
#             idx = int(mut[1:-1])-1
#             new_base = mut[-1]
#             muts[idx] = new_base
#         node_seq = edit_string(refseq, muts)
        
#         # Include each sequence once
#         if node_seq not in leaf2seq.values():
#             leaf2seq[leaf.id] = node_seq
    
#     out_path = f"{clade_path}/reconstructed_seqs.fasta"
#     leaf2seq["ancestral"] = refseq_orig
#     with open(out_path, 'w') as f:
#         for leaf_id, seq in leaf2seq.items():
#             f.write(f">{leaf_id}\n{seq.strip()}\n")

# def edit_string(original, mutations):
#     """
#     Given the original string and a list of mutations, returns the new edited string.
#     """
#     new_str = ""
#     prev_idx = 0
#     idx_list = list(mutations.keys())
#     idx_list.sort()
#     for idx in idx_list:
#         new_str += original[prev_idx:idx] + mutations[idx]
#         prev_idx = idx+1
#     if idx+1 < len(original):
#         new_str += original[prev_idx:len(original)+1]
    
#     # Check to make sure your new sequence is correct
#     assert len(new_str) == len(original)
#     for i in range(len(original)):
#         if i not in mutations.keys():
#             assert original[i] == new_str[i]
#         else:
#             new_str[i] == mutations[i]

#     return new_str
    


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
