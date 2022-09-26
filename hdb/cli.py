"""Command line interface."""

import sys

import click


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
@click.argument("settings_json", required=True, type=click.Path(exists=True))
@click.option("--make-paths-absolute", is_flag=True, help="Make paths absolute.")
def example(settings_json, make_paths_absolute):
    """Generate a file using one of our templates and the settings."""
    print("Invoked with", settings_json, make_paths_absolute)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
