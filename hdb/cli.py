"""Command line interface."""

import sys
import click
import hdb.summary as summary


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


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
