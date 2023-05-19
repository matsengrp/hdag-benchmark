import click
import random
import pickle
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import seaborn as sns
sns.set_theme()

import random
from math import exp
import math
from collections import Counter

from historydag import parsimony_utils
import historydag as hdag
from plot_utils import sliding_window_plot, get_results_full, bin_hist_plot

@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli():
    """
    Command line scripts to facilitate node support validation
    """
    pass

    
@click.command("coverage_trial_plot")
@click.option('--input', '-i', help='file path to input list as a pickle.')
@click.option('--out_dir', '-o', help='output directory to store figures/tables in.')
@click.option('--method', '-m', default='historydag')
@click.option('--clade_name', '-c')
@click.option('--window_proportion', '-w', default=0.20, help='the proportion of the data to use as window size')
def coverage_trial_plot(input, out_dir, clade_name, method, window_proportion=0.20):
    """
    Generates CA plot for a single trial.
    """

    try:
        with open(input, "rb") as f:
            results = pickle.load(f)
    except:
        print(f"\t...Skipping {input}")
        return

    # Remove leaf nodes
    res_no_leaves = []
    for el in results:
        if len(el[0]) > 1:
            res_no_leaves.append(el)
    results = res_no_leaves


    if method != "historydag":
        window_size = 20
    else:
        window_size = int(len(results) * window_proportion)

    out_path = out_dir + f"/support_quartiles_w={window_size}.png"
    x, y, min_sup, max_sup = sliding_window_plot(results, window_size=window_size, sup_range=True)

    # print("est\ttrue\tQ1\tQ3")
    # for num, (i, j, k, l) in enumerate(zip(x, y, min_sup, max_sup)):
    #     print(f"{num}\t{i:3f}\t{j:3f}\t{k:3f}\t{l:3f}")
    # print(window_size)

    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[0.75, 0.25])
    f.set_size_inches(7, 9)

    ax1.set_title(f"Clade {clade_name} Coverage Analysis")
    ax1.set_ylabel(f"Empirical Probability (window_size={window_size}/{len(results)})")

    ax1.plot(x, y, label="Support", color="red")
    ax1.scatter(min_sup, y, color="orange", alpha=0.1)
    ax1.scatter(max_sup, y, color="orange", alpha=0.1)
    # ax1.fill_betweenx(y, min_sup, max_sup, alpha=0.2, color="orange", label="Support Range")
    ax1.plot([0, 1], [0, 1], color="blue", label="Perfect", linestyle="dashed")
    ax1.legend()
    
    ax2.hist(x)
    ax2.set_xlabel("Estimated Support")
    ax2.set_yscale("log")
    f.savefig(out_path)
    f.clf()



@click.command("clade_results")
@click.option('--clade_dir', '-c', help='path to clade directory.')
@click.option('--out_path', '-o', help='output path to store figures/tables in.')
@click.option('--results_name', '-r', default="results.pkl", help='name of file (including path extension e.g. pkl).')
@click.option('--num_sim', '-n', default=1, help='number of simulations to average.')
@click.option('--trial_nums', '-t', help='the simulations to average over.')
@click.option('--method', '-m', default='historydag')
@click.option('--window_proportion', '-w', default=0.20, help='the proportion of the data to use as window size')
def clade_results(clade_dir, out_path, num_sim, method, window_proportion, results_name, trial_nums=None):
    """
    Generates CA plot by combining inferences across all trials for a given clade.

    Given the clade directory, performs coverage analysis across all simulations.
    """

    skip_list = []
    # print("Raw trial nums:", trial_nums)
    if trial_nums is not None:
        trial_num_ints = [int(n) for n in trial_nums.split(",")]
        # print("Trial nums:", trial_num_ints)
        skip_list = [i for i in range(1, num_sim+1) if i not in trial_num_ints]
        # print("Skipping:", skip_list)

    results_full = get_results_full(clade_dir, num_sim, method, results_name, skip_list)
    window_size = int(len(results_full) * window_proportion)

    print(f"\tgenerating window plot at {out_path}...")
    if method == "mrbayes":
        window_size = 100
    else:
        # window_size = int(len(results_full) * window_proportion)
        window_size=100

    # Sliding window plot
    x, y, pos_devs, neg_devs = sliding_window_plot(results_full, std_dev=True, window_size=window_size)

    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[0.75, 0.25])
    f.set_size_inches(7, 9)
    clade_name = clade_dir.split("/")[-1]
    ax1.set_title(f"Clade {clade_name} Coverage Analysis")

    ax1.set_ylabel(f"Empirical Probability (window_size={window_size}/{len(results_full)})")
    ax1.plot(x, y, color="orange", label="Support Regressor")
    ax1.fill_between(x, pos_devs, neg_devs, alpha=0.2, color="orange")
    ax1.plot([0, 1], [0, 1], color="blue", label="Perfect Regressor", linestyle="dashed")
    ax1.legend()
    
    ax2.hist(x)
    ax2.set_yscale("log")
    ax2.set_xlabel("Estimated Support")
    f.savefig(out_path)
    f.clf()

    # Binned histogram plot
    x, y, pos_devs, neg_devs = bin_hist_plot(results_full, bin_size=0.1)

    f, ax1 = plt.subplots(1, 1)
    f.set_size_inches(6.4, 4.8)
    ax1.set_title(f"Clade {clade_name} Binned Support")

    ax1.set_ylabel(f"Empirical Probability")
    ax1.plot(x, y, color="orange", label="Support Regressor")
    ax1.fill_between(x, pos_devs, neg_devs, alpha=0.2, color="orange")
    ax1.plot([0, 1], [0, 1], color="blue", label="Perfect Regressor", linestyle="dashed")
    ax1.legend()
    f.savefig(out_path[:-4]+"_histogram.png")
    f.clf()


@click.command("agg_pars_weights")
@click.option('--input', '-i', help='file path to input list as a pickle.')
@click.option('--out_dir', '-o', help='output directory to store figures/tables in.')
@click.option('--method', '-m', default='historydag')
@click.option('--clade_name', '-c')
@click.option('--window_proportion', '-w', default=0.20, help='the proportion of the data to use as window size')
def agg_pars_weights(input, out_dir, clade_name, method, window_proportion=0.20):
    """Given the pickled file, aggregates different inferences for various parsimony weights
    for computing support values
    
    E.g.
        python support_pipeline_scripts/cli.py agg_strats \
        -i /home/whowards/hdag-benchmark/data/A.2.2/1/results/historydag/strat_dict.pkl \
        -c A.2.2 \
        -o /home/whowards/hdag-benchmark/data/A.2.2/1/figures/historydag
    """

    try:
        with open(input, "rb") as f:
            strat_dict = pickle.load(f)
    except:
        print(f"\t...Skipping {input}")
        return

    # p_cmap = {
    #     0:#8dd3c7,
    #     1:#ffffb3,
    #     2:#bebada,
    #     3:#fb8072,
    #     4:#80b1d3,
    #     10:#fdb462,
    #     100:#b3de69,
    #     "inf":#fccde5
    # }

    # NOTE: Make sure the lengths match
    pweight = [0, 1, 2, 3, 4, 10, 1e2, 1e10, "inf"]
    colors = ['#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f', '#000000', '#33a02c']
    p_cmap = {p: c for p, c in zip(pweight, colors)}

    
    for strat, (p, results) in strat_dict.items():
        if p not in pweight:
            continue
        # Remove leaf nodes
        res_no_leaves = []
        for el in results:
            if len(el[0]) > 1:
                res_no_leaves.append(el)
        results = res_no_leaves


        window_size = int(len(results) * window_proportion)

        out_path = out_dir + f"/pars_weight_w={window_size}_complete.png"
        x, y, min_sup, max_sup = sliding_window_plot(results, window_size=window_size, sup_range=True)

        # f, ax1 = plt.subplots(1, 1, sharex=True)

        plt.plot(x, y, color=p_cmap[p], label=p)

        # TODO: Add these back in once you have a better idea of the correct pweight
        # plt.scatter(min_sup, y, color=p_cmap[p], alpha=0.1)
        # plt.scatter(max_sup, y, color=p_cmap[p], alpha=0.1)
        # plt.fill_betweenx(y, min_sup, max_sup, alpha=0.1, color=p_cmap[p])

    plt.title(f"Clade {clade_name} Coverage Analysis with Varying Parsimony Weights")
    plt.ylabel(f"Empirical Probability (window_size={window_size}/{len(results)})")
    plt.xlabel("Estimated Support")
    plt.plot([0, 1], [0, 1], color="blue", label="Perfect")
    # plt.legend()
    plt.savefig(out_path)
    plt.clf()


@click.command("bin_pars_weights")
@click.option('--input', '-i', help='file path to input list as a pickle.')
@click.option('--out_dir', '-o', help='output directory to store figures/tables in.')
@click.option('--method', '-m', default='historydag')
@click.option('--clade_name', '-c')
@click.option('--bin_size', '-b', default=0.05, help='the proportion of the data to use as window size')
def bin_pars_weights(input, out_dir, clade_name, method, bin_size):
    """Given the pickled file, aggregates different inferences for various parsimony weights
    for computing support values USING A HISTOGRAM BINNING STRATEGY (AS OPPOSED TO SLIDING WINDOW)
    
    E.g.:
    python support_pipeline_scripts/cli.py bin_pars_weights \
    -i /home/whowards/hdag-benchmark/data/A.2.2/1/results/historydag/strat_dict_pars_weight.pkl \
    -c A.2.2 \
    -o /home/whowards/hdag-benchmark/data/A.2.2/1/figures/historydag
    """

    try:
        with open(input, "rb") as f:
            strat_dict = pickle.load(f)
    except:
        print(f"\t...Skipping {input}")
        return

    pweight = [0, 1, 2, 3, 4, 10, 1e2, 1e10, "inf"]
    colors = ['#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f', '#000000', '#33a02c']
    p_cmap = {p: c for p, c in zip(pweight, colors)}

    
    for _, (p, results) in strat_dict.items():
        if p not in pweight:
            continue
        # Remove leaf nodes
        res_no_leaves = []
        for el in results:
            if len(el[0]) > 1:
                res_no_leaves.append(el)
        results = res_no_leaves

        # TODO: Change this to be for complete and complete_collapsed...
        out_path = out_dir + f"/pars_weight_binned_{bin_size}_collapse.png"
        x, y, std_pos, std_neg = bin_hist_plot(results, bin_size=bin_size)
        plt.plot(x, y, color=p_cmap[p], label=p)

        # TODO: Add these back in once you have a better idea of the correct pweight
        # plt.scatter(min_sup, y, color=p_cmap[p], alpha=0.1)
        # plt.scatter(max_sup, y, color=p_cmap[p], alpha=0.1)
        # plt.fill_betweenx(y, min_sup, max_sup, alpha=0.1, color=p_cmap[p])

    plt.title(f"Clade {clade_name} Coverage Analysis with Varying Parsimony Weights")
    plt.ylabel(f"Empirical Probability")
    plt.xlabel("Estimated Support")
    plt.plot([0, 1], [0, 1], color="blue", label="Perfect")
    # plt.legend()
    plt.savefig(out_path)
    plt.clf()


@click.command("agg_strats")
@click.option('--input', '-i', help='file path to input list as a pickle.')
@click.option('--out_dir', '-o', help='output directory to store figures/tables in.')
@click.option('--method', '-m', default='historydag')
@click.option('--clade_name', '-c')
@click.option('--window_proportion', '-w', default=0.20, help='the proportion of the data to use as window size')
def agg_strats(input, out_dir, clade_name, method, window_proportion=0.20):
    """Given the pickled file, aggregates different strategies for computing support values
    
    E.g.
        python support_pipeline_scripts/cli.py agg_strats \
        -i /home/whowards/hdag-benchmark/data/A.2.2/1/results/historydag/strat_dict.pkl \
        -c A.2.2 \
        -o /home/whowards/hdag-benchmark/data/A.2.2/1/figures/historydag

    """

    try:
        with open(input, "rb") as f:
            strat_dict = pickle.load(f)
    except:
        print(f"\t...Skipping {input}")
        return


    # TODO: Maybe we can use this type of plot to compare BEAST with MP-based support??

    mp_sups = []
    full_sups = []
    in_tree_labels = []

    mp_results = strat_dict["mp"][1]
    mp_dict = {node: (est_sup, in_tree) for node, est_sup, in_tree in mp_results}

    full_results = strat_dict["full"][1]
    for node, est_sup, in_tree in full_results:
        if node not in mp_dict:
            continue

        if len(node) == 1:
            print("Skipping leaf!")
            continue
        
        mp_sups.append(mp_dict[node][0])
        full_sups.append(est_sup)
        in_tree_labels.append(in_tree)

    x = [sup for sup, in_tree in zip(mp_sups, in_tree_labels) if in_tree]
    y = [sup for sup, in_tree in zip(full_sups, in_tree_labels) if in_tree]
    plt.scatter(x, y, label="In Tree", alpha=0.3)
    x = [sup for sup, in_tree in zip(mp_sups, in_tree_labels) if not in_tree]
    y = [sup for sup, in_tree in zip(full_sups, in_tree_labels) if not in_tree]
    plt.scatter(x, y, label="Not In Tree", alpha=0.3)
    plt.xlabel("MP-Trimmed Support")
    plt.ylabel("Full DAG Support")
    plt.title(f"Support Scatter")
    plt.plot([0, 1], [0, 1], color="green", linestyle="-", label="y = x", alpha=0.3)
    plt.legend()
    plt.savefig(out_dir + f"/mp_vs_full_support_comparison.png")
        

    
    for strat, (trimmed_weights, results) in strat_dict.items():
        # Remove leaf nodes
        res_no_leaves = []
        for el in results:
            if len(el[0]) > 1:
                res_no_leaves.append(el)
        results = res_no_leaves


        window_size = int(len(results) * window_proportion)

        out_path = out_dir + f"/support_quartiles_w={window_size}_strat={strat}.png"
        x, y, min_sup, max_sup = sliding_window_plot(results, window_size=window_size, sup_range=True)

        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[0.75, 0.25])
        f.set_size_inches(7, 10)

        ax1.set_title(f"Clade {clade_name} Coverage Analysis with \"{strat}\" Strategy")
        ax1.set_ylabel(f"Empirical Probability (window_size={window_size}/{len(results)})")

        ax1.plot(x, y, label="Support", color="red")
        ax1.scatter(min_sup, y, color="orange", alpha=0.1)
        ax1.scatter(max_sup, y, color="orange", alpha=0.1)
        # ax1.fill_betweenx(y, min_sup, max_sup, alpha=0.2, color="orange", label="Support Range")
        ax1.plot([0, 1], [0, 1], color="blue", label="Perfect")
        ax1.legend()
        
        ax2.hist(x)
        ax2.set_xlabel("Estimated Support")
        ax2.set_yscale("log")
        f.savefig(out_path)
        f.clf()

        f.set_size_inches(6, 8)
        plt.bar(list(trimmed_weights.keys()), list(trimmed_weights.values()))
        plt.xlabel("Parsimony")
        plt.yscale("log")
        plt.title(f"Parsiomny Distribution in DAG for \"{strat}\" Strategy")
        plt.savefig(out_dir + f"/parsiomny_distribution_for_{strat}.png")


@click.command("cumul_pars_weight")
@click.option('--input', '-i', help="path to input MADAG")
@click.option('--out_dir', '-o', help="directory path to output plot to")
@click.option('--parsimony_weight', '-p', default=0.01, help="the coefficient to multiple parsiomny by in the negative exponential")
def cumul_pars_weight(input, out_dir, parsimony_weight):
    """
    Plots the cummulative distribution of tree probability as a function of parsiomny score
    for various parsimony weightings.

    Roughly, this gives us a sense of the probability of sampling a tree with a given parsimony score or lower.
    """

    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(input)
    pars_distribution = dag.weight_count()

    p_scores = list(pars_distribution.keys())
    p_scores.sort()
    cumm_probs = {}
    print("\t Number of Parsimony scores:", len(p_scores))
    for k in [-2, -1, -0.5, 0, 0.1, 0.5, 1]:
        print("\t", k)
        # dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(input) # TODO: Try getting a fresh copy of the dag.. Doesn't help
        total_tree_score = dag.probability_annotate(lambda n1, n2: -k * parsimony_utils.hamming_cg_edge_weight(n1, n2), log_probabilities=True)
        cumm_probs[k] = []
        for p in p_scores:
            val = pars_distribution[p] * math.e ** (-k * p) / exp(total_tree_score)
            if len(cumm_probs[k]) < 1:
                cumm_probs[k].append(val)
            else:
                cumm_probs[k].append(val + cumm_probs[k][-1])

    for k, cumm_prob in cumm_probs.items():
        plt.plot(p_scores, cumm_probs[k], label=k)

    plt.legend()
    plt.xlabel("Parsimony Score")
    plt.ylabel("Cummulative Probability Weight")
    plt.title(f"Cummulative Probability with Parsimony Weight {parsimony_weight}")
    out_path = out_dir + f"/cumm_prob_k={parsimony_weight}.png"
    plt.savefig(out_path)


@click.command("clade_results_random_scaling")
@click.option('--clade_dir', '-c', help='path to clade directory.')
@click.option('--out_path', '-o', help='output path to store figures/tables in.')
@click.option('--results_name', '-r', default="results.pkl", help='name of file (including path extension e.g. pkl).')
@click.option('--num_sim', '-n', default=1, help='number of simulations to average.')
@click.option('--method', '-m', default='historydag')
@click.option('--sample_size', '-s', default=0.5, help='proportion of results to randomly scale.')
@click.option('--window_proportion', '-w', default=0.20, help='the proportion of the data to use as window size')
def clade_results_random_scaling(clade_dir, out_path, num_sim, method, window_proportion, results_name, sample_size):
    """Given the clade directory, performs coverage analysis across all simulations"""

    # scale_range = (0.75, 0.95)
    scale_range = (0.5, 1)

    # Randomly scale support values and re-sort
    results_full = get_results_full(clade_dir, num_sim, method, results_name)
    print("length of results", len(results_full))
    results_full_rand = []
    idxs = random.sample(range(len(results_full)), int(len(results_full) * sample_size))
    for i, (node, est_sup, in_tree) in enumerate(results_full):
        if i in idxs:
            est_sup = est_sup * random.uniform(scale_range[0], scale_range[1])
        results_full_rand.append((node, est_sup, in_tree))
    results_full = results_full_rand
    random.shuffle(results_full)
    results_full.sort(key=lambda el: el[1])

    if method == "mrbayes":
        window_size = 500
    else:
        window_size = int(len(results_full) * window_proportion)

    print(f"\tgenerating window plot at {out_path}...")
    x, y, pos_devs, neg_devs = sliding_window_plot(results_full, std_dev=True, window_size=window_size)

    # Plotting code
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[0.75, 0.25])
    f.set_size_inches(7, 9)
    clade_name = clade_dir.split("/")[-1]
    ax1.set_title(f"Clade {clade_name} Coverage Analysis")
    ax1.set_ylabel(f"Empirical Probability (window_size={window_size}/{len(results_full)})")
    ax1.plot(x, y, color="orange", label="Support Regressor")
    ax1.fill_between(x, pos_devs, neg_devs, alpha=0.2, color="orange")
    ax1.plot([0, 1], [0, 1], color="blue", label="Perfect Regressor", linestyle="dashed")
    ax1.legend()
    ax2.hist(x)
    ax2.set_yscale("log")
    ax2.set_xlabel("Estimated Support")
    f.savefig(out_path)
    f.clf()


# TODO: Rename this. This plots all parsiomny coefficient coverage analyses on the same plot with
#       a sliding window.
@click.command("pars_weight_clade_results")
@click.option('--clade_dir', '-c', help='path to clade directory.')
@click.option('--out_dir', '-o', help='output directory to store figures/tables in.')
@click.option('--num_sim', '-n', default=1, help='number of simulations to average.')
@click.option('--method', '-m', default='historydag')
@click.option('--bin_size', '-b', default=0.05, help='the proportion of the data to use as window size')
def pars_weight_clade_results(clade_dir, out_dir, num_sim, method, bin_size):
    """Given the clade directory with a results dictionary, performs coverage analysis across all simulations"""

    clade_name = clade_dir.split("/")[-1]

    result_dict = {}
    for trial in range(1, num_sim+1):
        # Assumes that `path/to/clade/trial/results/method/strat_dict_node_weight.pkl stores``
        #   list of nodes their supports and whether they're in the true tree or not
        result_path = clade_dir + f"/{trial}/results/{method}/strat_dict_node_weight.pkl"
        
        try:
            with open(result_path, "rb") as f:
                results = pickle.load(f)
                result_dict[trial] = results
        except:
            print(f"\tSkipping {result_path}")
            continue
    
    if len(result_dict) == 0:
        print("\n==> No results to print :(\n")
        return

    results_full = {p: [] for p in result_dict[1].keys() if isinstance(p, str) or p <= 4}
    for trial, strat_dict in result_dict.items():
        # NOTE: First element of results should be p and second element is the results list
        for p, results in strat_dict.items():
            if p not in results_full:
                continue

            results = results[1]
            
            # NOTE: Removing all leaves here
            print(f"p={p}, trial={trial}")
            with_leaves = len(results)
            leaf_in_tree = [int(result[2]) for result in results if len(result[0]) <= 1]
            leaf_est_sup = [result[1] for result in results if len(result[0]) <= 1]
            results = [result for result in results if len(result[0]) > 1 and result[1] > 0]    # NOTE: Only including nodes in DAG!
            without_leaves = len(results)
            if with_leaves != without_leaves:
                print(f"==> Removed {with_leaves - without_leaves} nodes \
                    avg in_tree = {sum(leaf_in_tree) / len(leaf_in_tree)} \
                    avg est_sup = {sum(leaf_est_sup) / len(leaf_est_sup)}")
            ##

            results_full[p].extend(results)

    # print(f"\tsorting {len(results_full)} results...")
    for _, results in results_full.items():
        random.shuffle(results)
        results.sort(key=lambda el: el[1])

    bin_size = 0.05
    plt.plot([0, 1], [0, 1], color="blue", label="Perfect", linestyle='dashed')

    pweight = list(results_full.keys())
    for val in ["inf", "bif"]:
        if val in pweight:
            pweight.remove(val)

    base = 10
    print(pweight)
    colors = list(plt.cm.autumn((np.power(base, np.linspace(0, 1, int(max(pweight) * 10) + 1)) - 1) / (base - 1)))
    colors_neg = list(plt.cm.winter((np.power(base, np.linspace(0, 1, int(-1 * min(pweight) * 10) + 1)) - 1) / (base - 1)))
    colors.reverse()
    # colors_neg.reverse()

    for p, results in results_full.items():
        if not isinstance(p, str) and p > 10:
            continue
        # Remove leaf nodes
        res_no_leaves = []
        for el in results:
            if len(el[0]) > 1:
                res_no_leaves.append(el)
        results = res_no_leaves

        out_path = out_dir + f"/single_line_pars_weighted_w={200}.png"
        x, y, std_pos, std_neg = sliding_window_plot(results, std_dev=True, window_size=200)

        if not isinstance(p, str) and p >= 0:
            c = colors[int(p*10)]
        elif not isinstance(p, str) and p < 0:
            c = colors_neg[int(-1*p*10)]
        elif p == "inf":
            c = "green"
        else:
            c = "purple"
        
        plt.plot(x, y, color=c, label=p)

        # TODO: Add these back in once you have a better idea of the correct pweight
        # plt.scatter(min_sup, y, color=p_cmap[p], alpha=0.1)
        # plt.scatter(max_sup, y, color=p_cmap[p], alpha=0.1)
        # plt.fill_betweenx(y, min_sup, max_sup, alpha=0.1, color=p_cmap[p])

    plt.title(f"Clade {clade_name} Coverage Analysis with Varying Parsimony Weights")
    plt.ylabel(f"Empirical Probability")
    plt.xlabel("Estimated Support")
    plt.legend()
    plt.savefig(out_path)
    plt.clf()

    # Plot stuff about the MP results
    # for p in results_full.keys():
    #     results = results_full[p]

    #     num_nodes_in_tree_but_not_dag = sum([result[2] and result[1] == 0 for result in results])
    #     print("num_nodes_in_tree_but_not_dag =", num_nodes_in_tree_but_not_dag)

    #     # Support distribution
    #     in_tree_sup = [result[1] for result in results if result[2]]
    #     out_tree_sup = [result[1] for result in results if not result[2]]


    #     plt.hist((in_tree_sup, out_tree_sup), label=("In True Tree","Out True Tree"), color=("Orange", "Blue"))
    #     plt.ylabel("Count")
    #     plt.xlabel("Estimated Support")
    #     plt.legend()
    #     plt.savefig(out_dir + f"/{p}_support_histogram.png")
    #     plt.clf()

    #     # Scatter support vs clade size
    #     in_tree_size = [len(result[0]) for result in results if result[2]]
    #     out_tree_size = [len(result[0]) for result in results if not result[2]]

    #     plt.scatter(out_tree_size, out_tree_sup, alpha=0.15, c="Blue", label="Out True Tree")
    #     plt.scatter(in_tree_size, in_tree_sup, alpha=0.3, c="Orange", label="In True Tree")
    #     plt.ylabel("Estimated Support")
    #     plt.xlabel("Size")
    #     plt.legend()
    #     plt.savefig(out_dir + f"/{p}_support_size_scatter.png")
    #     plt.clf()

cli.add_command(coverage_trial_plot)
cli.add_command(clade_results)

    

if __name__ == '__main__':
    cli()
