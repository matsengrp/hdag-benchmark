"""Command line interface."""

import sys
import click
import hdb.summary as summary
import hdb.collapse_tree as ct


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
        intree = ete3.Tree(fh.read())
    outtree = ct.collapse_tree(intree, load_fasta(input_fasta))
    outtree.write(features=["mutations"], format_root_node=True, outfile=output_newick)

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
