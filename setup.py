"""Our setup script."""

import glob
from setuptools import setup

setup(
    name="hdb",
    description="Utilities for history DAG benchmarking.",
    packages=["hdb"],
    package_data={"hdb": ["data/*"]},
    scripts=glob.glob("hdb/scripts/*.sh"),
    entry_points={"console_scripts": ["hdb=hdb.cli:safe_cli"]},
    install_requires=[
        "PyQt5",
        "amas",
        "biopython",
        "click",
        "ete3",
        "historydag",
        "matplotlib",
        "pandas",
        "phastSim",
        "seqmagick"
    ],
)
