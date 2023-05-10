import click
import historydag as hdag
import ete3 as ete
import random
import pickle
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import dendropy
import json

import random
import subprocess
from math import exp
import math
import time
from collections import Counter

# TODO: Uncomment this later
from historydag import parsimony_utils
from historydag.parsimony import parsimony_score, sankoff_upward
from historydag.utils import count_labeled_binary_topologies

import seaborn as sns
sns.set_theme()
# plt.rcParams['text.usetex'] = True

@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Command line scripts to facilitate node support validation
    """
    pass

def support_adjustment():
    pass


cli.add_command(support_adjustment)

if __name__ == '__main__':
    cli()