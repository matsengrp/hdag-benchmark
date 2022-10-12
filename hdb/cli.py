"""Command line interface."""

import sys
import click
import hdb.summary as summary
import hdb.collapse_tree as ct
import hdb.aggregate_dnapars_trees as agg
import ete3
import pickle


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


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
