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


@cli.command()
@click.option("-i", "--input-path", help="Newick tree input path.")
@click.option("-o", "--output-path", help="Output path.")
@click.option("-s", "--resolve-seed", default=1)
@click.option("-b", "--branch-len-model", default="num-muts")
def resolve_multifurcations(input_path, output_path, resolve_seed, branch_len_model):
    """
    Given an ete tree, resolves all polytomies by creating a
    uniformly random bifurcating tree that is consistent with
    the multifurcating one.
    """
    tree = ete3.Tree(input_path, format=1)
    resolve_polytomy(tree, resolve_seed, branch_len_model)
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


@cli.command()
@click.option("-c", "--clade-path", help="Path to clade directory.")
@click.option("-n", "--num-trials", default=25, help="Number of trials to return.")
def find_diversity(clade_path, num_trials):
    trial_list = []
    for trial in range(1, 101):
        log_path = clade_path + f"/{trial}/results/historydag/opt_info/optimization_log_complete/logfile.csv"
        try:
            with open(log_path, "r") as f:
                num_trees = f.readlines()[-1].split("\t")[-3]
                trial_list.append((trial, num_trees))
        except:
            print(f"\tSkipping {log_path}")
            continue
    
    trial_list.sort(key=lambda el: el[1])
    trial_list = trial_list[-num_trials:]
    with open("pars_div_trials.txt", "w") as f:
        for trial, _ in trial_list:
            f.write(f"{trial}\n")
    


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
