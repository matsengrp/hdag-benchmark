import click
import random
import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import json
from scipy import stats
from scipy.stats import linregress, t
import pandas as pd

import random
from math import exp
import math
from collections import Counter

# from historydag import parsimony_utils
import historydag as hdag
from plot_utils import sliding_window_plot, get_results_general, bin_hist_plot, plot_scatter_with_line

import matplotlib as mpl
import matplotlib.font_manager
import seaborn as sns
sns.set_theme()
sns.set_style("white")

# mpl.rcParams.update({
#     'font.size': 16,
#     'axes.titlesize': 17,
#     'axes.labelsize': 16,
#     'xtick.labelsize': 13,
#     'ytick.labelsize': 13,
#     'font.family': 'sans-serif',
#     'font.weight': 600,
#     'axes.labelweight': 600,
#     'axes.titleweight': 600,
#     'figure.autolayout': True
#     })


colors = [
    '#66c2a5',
    '#fc8d62',
    '#8da0cb',
    '#e78ac3'
]


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

    cleaned_results = []
    for el in results:
        if len(el[0]) > 1: # Remove leaf nodes
            if el[1] > 0.01: # Remove low support nodes for efficiency
                cleaned_results.append(el)
    results = cleaned_results

    # print("DEBUG: Here are the supports your workign with")
    # print([result[1] for result in results])


    if method == "mrbayes":
        window_size = 100
    else:
        window_size = min(int(len(results) * window_proportion), 100)

    out_path = out_dir + f"/support_quartiles_w={window_size}.png"
    x, y, min_sup, max_sup = sliding_window_plot(results, window_size=window_size, sup_range=True)

    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[0.75, 0.25])
    f.set_size_inches(7, 9)

    ax1.set_title(f"Clade {clade_name} Coverage Analysis")
    ax1.set_ylabel(f"Empirical Probability (window_size={window_size}/{len(results)})")

    ax1.plot(x, y, label="Support", color="red")
    ax1.scatter(min_sup, y, color="orange", alpha=0.1)
    ax1.scatter(max_sup, y, color="orange", alpha=0.1)
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
@click.option('--trial_nums', '-t', default=None, help='the simulations to average over.')
@click.option('--method', '-m', default='historydag')
@click.option('--window_proportion', '-w', default=0.20, help='the proportion of the data to use as window size')
def clade_results(clade_dir, out_path, num_sim, method, window_proportion, results_name, trial_nums):
    """
    Generates CA plot by combining inferences across all trials for a given clade.

    Given the clade directory, performs coverage analysis across all simulations.
    """

    skip_list = []
    if trial_nums is not None:
        trial_num_ints = [int(n) for n in trial_nums.split(",")]
        skip_list = [i for i in range(1, num_sim+1) if i not in trial_num_ints]

    print("Gathering results...")
    results_full = get_results_general(clade_dir, num_sim, method, results_name, skip_list, support_removal_threshold=0.01)
    print([(el[1], el[2]) for el in results_full[-100:]])
    window_size = int(len(results_full) * window_proportion)

    print(f"\tgenerating window plot at {out_path}...")
    if method == "mrbayes":
        window_size = 100
    else:
        # window_size = int(len(results_full) * window_proportion)
        window_size = 100

    # Sliding window plot
    x, y, pos_devs, neg_devs = sliding_window_plot(results_full, std_dev=True, window_size=window_size)

    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[0.75, 0.25])
    f.set_size_inches(7, 9)
    clade_name = clade_dir.split("/")[-1] # TODO: This is a terrible way to get the clade name. Change it.
    ax1.set_title(f"Clade {clade_name} Coverage Analysis")

    ax1.set_ylabel(f"Empirical Probability (window_size={window_size}/{len(results_full)})")
    ax1.plot(x, y, color="orange", label="Support Regressor")
    ax1.fill_between(x, pos_devs, neg_devs, alpha=0.2, color="orange")
    ax1.plot([0, 1], [0, 1], color="blue", label="Perfect Regressor", linestyle="dashed")
    ax1.legend()
    
    ax2.hist(x)
    ax2.set_yscale("log")
    ax2.set_xlabel("Estimated Support")
    f.savefig(out_path + f"_{window_size}_.png")
    f.clf()

    # # Binned histogram plot
    # x, y, pos_devs, neg_devs = bin_hist_plot(results_full, bin_size=0.1)

    # f, ax1 = plt.subplots(1, 1)
    # f.set_size_inches(6.4, 4.8)
    # ax1.set_title(f"Clade {clade_name} Binned Support")

    # ax1.set_ylabel(f"Empirical Probability")
    # ax1.plot(x, y, color="orange", label="Support Regressor")
    # ax1.fill_between(x, pos_devs, neg_devs, alpha=0.2, color="orange")
    # ax1.plot([0, 1], [0, 1], color="blue", label="Perfect Regressor", linestyle="dashed")
    # ax1.legend()
    # f.savefig(out_path[:-4]+"_histogram.png")
    # f.clf()

@click.command("compare_trial_results")
@click.option('--clade_dir', '-c', help='path to clade directory.')
@click.option('--out_path', '-o', help='output path to store figures/tables in.')
@click.option('--results_name', '-r', default="results.pkl", help='name of file (including path extension e.g. pkl).')
@click.option('--trial_nums', '-t', help='the simulations to average over.')
@click.option('--method1', '-m1', default='historydag')
@click.option('--method2', '-m2', default='mrbayes')
def compare_trial_results(clade_dir, out_path, method1, method2, results_name, trial_nums=None):
    """
    Given two methods and the name of their results files, plots a scatter plot of their
    estimated supports.
    
    Uses only the nodes in method1 #TODO: Might wanna change this

    E.g.,
        python support_pipeline/plotting.py compare_trial_results \
        -c /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.108 \
        -o /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.108/figures/clade_comparison_trials \
        -t 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25
    
    """
    num_sim = 25
    skip_list = []
    # print("Raw trial nums:", trial_nums)
    if trial_nums is not None:
        trial_num_ints = [int(n) for n in trial_nums.split(",")]
        skip_list = [i for i in range(1, num_sim+1) if i not in trial_num_ints]

    for trial in range(1, num_sim+1):
        result_path1 = clade_dir + f"/{trial}/results/{method1}/{results_name}"
        result_path2 = clade_dir + f"/{trial}/results/{method2}/{results_name}"

        # Skip any trials you don't want (e.g., inference on them hasn't finished for this run)
        if trial in skip_list:
            continue
        
        try:
            # Convert lists to map of clade -> (est_support, in_tree)
            with open(result_path1, "rb") as f:
                results1 = pickle.load(f)
                results1 = {result[0]: (result[1], result[2]) for result in results1 if len(result[0]) > 1}
            with open(result_path2, "rb") as f:
                results2 = pickle.load(f)
                results2 = {result[0]: (result[1], result[2]) for result in results2 if len(result[0]) > 1}

            # Create list (node, est1, est2, in_tree)
            paired_results = []
            for node1 in results1.keys():
                in_tree = results1[node1][1]
                est1 = results1[node1][0]
                if node1 in results2:
                    est2 = results2[node1][0]
                else:
                    est2 = 0
                paired_results.append((node1, est1, est2, in_tree))
                
            print(trial)

            a=-1
            x_in = []
            y_in = []
            x_in.extend([est1 for node, est1, est2, in_tree in paired_results if in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a])
            y_in.extend([est2 for node, est1, est2, in_tree in paired_results if in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a])
            x_out = []
            y_out = []
            x_out.extend([est1 for node, est1, est2, in_tree in paired_results if not in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a])
            y_out.extend([est2 for node, est1, est2, in_tree in paired_results if not in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a])

            x = x_in.copy()
            x.extend(x_out)
            y = y_in.copy()
            y.extend(y_out)

            data = np.array([x, y]).T
            df = pd.DataFrame(data, columns = [method1, method2])
            plot_scatter_with_line(df, method1, method2, out_path+f'/{trial}_{method1}_{method2}.png')
            
            correlation = stats.pearsonr(y, x)
            print("\tCorrelation:", correlation)
            
            plt.clf()
            plt.scatter(x_in, y_in, alpha = 0.5, label="In tree")
            plt.scatter(x_out, y_out, alpha = 0.3, label="Not in tree")
            plt.xlabel(f"{method1} estimated supports")
            plt.ylabel(f"{method2} estimated supports")
            plt.legend()
            plt.savefig(out_path + f"/{trial}_{method1}_{method2}_scatter.png")

        except:
            print(f"---- Skipping {clade_dir} {trial} ----")
            continue
    


@click.command("compare_clade_results")
@click.option('--clade_dir', '-c', help='path to clade directory.')
@click.option('--out_path', '-o', help='output path to store figures/tables in.')
@click.option('--results_name', '-r', default="results.pkl", help='name of file (including path extension e.g. pkl).')
@click.option('--trial_nums', '-t', help='the simulations to average over.')
@click.option('--method1', '-m1', default='historydag')
@click.option('--method2', '-m2', default='mrbayes')
def compare_clade_results(clade_dir, out_path, method1, method2, results_name, trial_nums=None):
    """
    Given two methods and the name of their results files, plots a scatter plot of their
    estimated supports.
    
    Uses only the nodes in method1 #TODO: Might wanna change this

    E.g.,
        python support_pipeline/plotting.py compare_clade_results \
        -c /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.108 \
        -o /fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.108/figures \
        -t 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25
    
    """
    num_sim = 25
    skip_list = []
    # print("Raw trial nums:", trial_nums)
    if trial_nums is not None:
        trial_num_ints = [int(n) for n in trial_nums.split(",")]
        skip_list = [i for i in range(1, num_sim+1) if i not in trial_num_ints]

    result_dict = {}
    for trial in range(1, num_sim+1):
        # TODO: Figure out how to not include trials in a more natural way
        result_path1 = clade_dir + f"/{trial}/results/{method1}/{results_name}"
        result_path2 = clade_dir + f"/{trial}/results/{method2}/{results_name}"

        # Skip any trials you don't want (e.g., inference on them hasn't finished for this run)
        if trial in skip_list:
            continue

        
        try:
            # Convert lists to map of clade -> (est_support, in_tree)
            with open(result_path1, "rb") as f:
                results1 = pickle.load(f)
                results1 = {result[0]: (result[1], result[2]) for result in results1 if len(result[0]) > 1}
            with open(result_path2, "rb") as f:
                results2 = pickle.load(f)
                results2 = {result[0]: (result[1], result[2]) for result in results2 if len(result[0]) > 1}

            # Create list (node, est1, est2, in_tree)
            paired_results = []
            for node1 in results1.keys():
                in_tree = results1[node1][1]
                est1 = results1[node1][0]
                if node1 in results2:
                    est2 = results2[node1][0]
                else:
                    est2 = 0
                paired_results.append((node1, est1, est2, in_tree))
                
            result_dict[trial] = paired_results
            print(trial)
        except:
            print(f"\t--- Skipping {clade_dir} {trial} ---")
            continue
    
    if len(result_dict) == 0:
        print("\n==>No results to print :(\n")
        return
    
    with open(f"{out_path}_comparsion_{method1}_{method2}.pkl", "wb") as f:
        pickle.dump(result_dict, f)

    with open(f"{out_path}_comparsion_{method1}_{method2}.pkl", "rb") as f:
        result_dict = pickle.load(f)

    a=0.025
    x_in = []
    y_in = []
    for trial, results in result_dict.items():
        x_in.extend([est1 for node, est1, est2, in_tree in results if in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a])
        y_in.extend([est2 for node, est1, est2, in_tree in results if in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a])
    x_out = []
    y_out = []
    for trial, results in result_dict.items():
        x_out.extend([est1 for node, est1, est2, in_tree in results if not in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a])
        y_out.extend([est2 for node, est1, est2, in_tree in results if not in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a])

    x = x_in.copy()
    x.extend(x_out)
    y = y_in.copy()
    y.extend(y_out)

    data = np.array([x, y]).T
    df = pd.DataFrame(data, columns = [method1, method2])
    plot_scatter_with_line(df, method1, method2, out_path+f'/clade_comparison_{method1}_{method2}.png')
    
    correlation = stats.pearsonr(y, x)
    print("Correlation:", correlation)
    
    
    # plt.scatter(x_in, y_in, alpha = 0.3, label="In tree")
    # plt.scatter(x_out, y_out, alpha = 0.3, label="Not in tree")
    # plt.xlabel(f"{method1} estimated supports")
    # plt.ylabel(f"{method2} estimated supports")
    # plt.legend()
    # plt.savefig(out_path + f"/clade_comparison_{method1}_{method2}.png")


@click.command("aggregate_results")
@click.option('--base_dir', '-c', help='path to base directory where all clades reside.')
@click.option('--out_dir', '-o', help='output path to store figures/tables in.')
@click.option('--results_name', '-r', default="results.pkl", help='name of file (including path extension e.g. pkl).')
@click.option('--method', '-m', default='historydag,')
@click.option('--support_removal_threshold', '-t', default=0.01)
def aggregate_results(base_dir, out_dir, method, results_name, support_removal_threshold):
    """
    Generates CA plot by combining inferences across all clades for a given directory.

    E.g.
python support_pipeline/plotting.py aggregate_results \
-c /fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models \
-o /fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models/figures \
-r results.pkl \
-t 0.01 \
-m historydag,diff_historydag
    """
    method_str = method
    if ',' not in method:
        methods = [method]
    else:
        methods = method.split(",")
        if len(methods[-1]) == 0:
            methods = methods[:-1]

    method2color = {m: colors[i] for i, m in enumerate(methods)}

    clade_names = ['AY.34.2', 'AZ.3', 'AY.108', 'A.2.5', 'AY.87', 'AY.74', 'B.1.1.10', 'P.1.7']
    # clade_names = [f"{c}_" for c in clade_names]  # TODO: Change names so you don't have to do this for the real data

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[0.75, 0.25])
    fig.set_size_inches(7, 9)
    ax1.set_title(f"Coverage Analysis")
    ax1.set_ylabel(f"Percentage Correct")
    hist = []

    trials = range(25)
    
    print(f"generating window plot at {out_dir}...")
    for method in methods:
        print(f"Collecting results for {method}...")
        results_full = []
        for clade in clade_names:
            for trial in trials:
                fp = f"{base_dir}/{clade}/gamma_10_hmut_50/{trial}/results/{method}/{results_name}"
                try:
                    with open(fp, "rb") as f:
                        results = pickle.load(f)
                        cleaned_results = [result for result in results if len(result[0]) > 1 and result[1] > support_removal_threshold]
                        
                        # TODO: Adjust minimum cleaned results requirements...
                        # if len(cleaned_results) > 90:

                        print(clade, "\t with", len(cleaned_results))
                        results_full.extend(cleaned_results)
                except:
                    continue
        
        print(f"Building sliding window plot...")
        results_full.sort(key=lambda result: result[1])


        if method == "mrbayes":
            window_size = 100
        else:
            # window_size = int(len(results_full) * window_proportion)
            window_size = 100

        x, y, pos_devs, neg_devs = sliding_window_plot(results_full, std_dev=True, window_size=window_size)
        ax1.plot(x, y, color=method2color[method], label=f"{method}")
        ax1.fill_between(x, pos_devs, neg_devs, alpha=0.1, color=method2color[method])

        hist.append(x)

    ax2.hist(hist, color=list(method2color.values()))

    ax1.plot([0, 1], [0, 1], color="blue", alpha=0.5, label="ideal", linestyle="dashed")
    ax1.legend()
    
    ax2.set_yscale("log")
    ax2.set_xlabel("Estimated Support")
    fig.savefig(out_dir + f"/clade_comparison_{method_str}")
    fig.clf()


@click.command("compare_results")
@click.option('--base_dir', '-c', help='path to clade directory.')
@click.option('--out_path', '-o', help='output path to store figures/tables in.')
@click.option('--results_name', '-r', default="results.pkl", help='name of file (including path extension e.g. pkl).')
@click.option('--method1', '-m1', default='historydag')
@click.option('--method2', '-m2', default='mrbayes')
@click.option('--support_removal_threshold', '-t', default=0.01)
# @click.option('--rerun_results', , default=0.01)
def compare_results(base_dir, out_path, method1, method2, results_name, support_removal_threshold):
    """
    Given two methods and the name of their results files, plots a scatter plot of the
    estimated supports they give for each node in method1. Will aggregate the plot over all clades.
    
    Uses only the nodes in method1 #TODO: Might wanna change this

    E.g.,
        python support_pipeline/plotting.py compare_results \
        -c /fh/fast/matsen_e/whowards/hdag-benchmark/data \
        -o /fh/fast/matsen_e/whowards/hdag-benchmark/data/figures \
        -r results_adj.pkl \
        -m1 historydag \
        -m2 mrbayes
    """

    # clade_names = ['AY.34.2', 'AZ.3', 'AY.108', 'B.1.258.3', 'BA.1.5', 'B.1.1.432', 'AY.32', 'P.1.7']
    clade_names = ['AY.34.2', 'AZ.3', 'AY.108', 'B.1.1.432', 'P.1.7', 'AY.32']   # NOTE: These are the clades that are large enough
    clade_names = [f"{c}_" for c in clade_names]

    rerun_results=False

    if rerun_results:
        result_dict = {}
        for clade in clade_names:
            result_path1 = base_dir + f"/{clade}/results/{method1}/{results_name}"
            result_path2 = base_dir + f"/{clade}/results/{method2}/results.pkl"

            # print(result_path1)
            # print(result_path2)

            try:
                # Convert lists to map of clade -> (est_support, in_tree)
                with open(result_path1, "rb") as f:
                    results1 = pickle.load(f)
                    results1 = {result[0]: (result[1], result[2]) for result in results1 if len(result[0]) > 1 and result[1] > support_removal_threshold}
                with open(result_path2, "rb") as f:
                    results2 = pickle.load(f)
                    results2 = {result[0]: (result[1], result[2]) for result in results2 if len(result[0]) > 1 and result[1] > support_removal_threshold}

                # Create list (node, est1, est2, in_tree)
                paired_results = []
                for node1 in results1.keys():
                    in_tree = results1[node1][1]
                    est1 = results1[node1][0]
                    if node1 in results2:
                        est2 = results2[node1][0]
                    else:
                        est2 = 0
                    paired_results.append((node1, est1, est2, in_tree))

                if len(paired_results) > 0:
                    result_dict[clade] = paired_results
                    print(clade)
                else:
                    raise(Warning(f"Clade only has nonzero support for {len(paired_results)} nodes."))
            except:
                print(f"\t--- Skipping {clade} ---")
                continue
        
        if len(result_dict) == 0:
            print("\n==>No results to print :(\n")
            return
        
        with open(f"{out_path}/comparsion_{method1}_{method2}.pkl", "wb") as f:
            pickle.dump(result_dict, f)

    with open(f"{out_path}/comparsion_{method1}_{method2}.pkl", "rb") as f:
        result_dict = pickle.load(f)

    a=-1
    x_in = []
    y_in = []
    x_out = []
    y_out = []
    for clade, results in result_dict.items():
        # if clade != 'P.1.7_':
        #     continue
        clade_results1 = [est1 for node, est1, est2, in_tree in results if in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a]
        clade_results2 = [est2 for node, est1, est2, in_tree in results if in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a]
        x_in.extend(clade_results1)
        y_in.extend(clade_results2)
        print(f"{clade}\t supports {len(clade_results1)} nodes")
        x_out.extend([est1 for node, est1, est2, in_tree in results if not in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a])
        y_out.extend([est2 for node, est1, est2, in_tree in results if not in_tree and est1 > a and est2 > a and est1 < 1-a and est2 < 1-a])

    x = x_in.copy()
    x.extend(x_out)
    y = y_in.copy()
    y.extend(y_out)

    data = np.array([x, y]).T
    df = pd.DataFrame(data, columns = [method1, method2])
    plot_scatter_with_line(df, method1, method2, out_path+f'/scatter_{method1}_{method2}.png')
    
    correlation = stats.pearsonr(y, x)
    print("Correlation:", correlation)


#################################################
### Hypermutation Experiments
#################################################

@click.command("plot_hmut")
@click.option("-d", "--data-path", help="Path to data directory.")
def plot_hmut(data_path):
    """
    Gathers statistics about hypermutation simulations and creates figures summarizing the results.
    """

    clade_names = [
        "AY.34.2",
        "AY.108",
        "AZ.3",
        "P.1.7"
    ]

    hmut_params = [
        (0.01, 200),
        (0.001, 200),
        (0.0001, 200),
        (0.01, 20),
        (0.001, 20),
        (0.0001, 20),
        (0.01, 100),
        (0.001, 100),
        (0.0001, 100),
        (0.01, 10),
        (0.001, 10),
        (0.0001, 10)
    ]

    trials = [1, 2, 3]

    param_trial_list = [(el1, el2, el3) for el1 in clade_names for el2 in hmut_params for el3 in trials]

    trial_list = []

    # Iterate over every combination of hmut param and trial
    for clade_name, (prob, rate), trial in param_trial_list:
        base_path = f"{data_path}/{clade_name}/{prob}_{rate}/{trial}"

        log_path = base_path + "/results/historydag/opt_info/optimization_log_complete/logfile.csv"
        sim_info_path = base_path + f"/simulation/tree_stats.json"
        results_path = base_path + f"/results/historydag/results.pkl"
        support_log_path = base_path + f"/results/inference_historydag.log"
        real_data_results_path = data_path + f"/../real_data/{clade_name}/results/historydag/results.pkl"
        real_data_hdag_log_path = data_path + f"/../real_data/{clade_name}/results/historydag/opt_info/optimization_log_complete/logfile.csv"
        real_data_fasta_path = data_path + f"/../real_data/{clade_name}/reconstructed_seqs.fasta"

        try:
            with open(real_data_hdag_log_path, "r") as f:
                real_mp_score = int(f.readlines()[-1].split("\t")[4])

            with open(real_data_fasta_path, "r") as f:
                real_num_leaves = int(len(f.readlines())/2)

            with open(log_path, "r") as f:
                mp_score = int(f.readlines()[-1].split("\t")[4])

            with open(results_path, "rb") as f:
                results = pickle.load(f)
                results_in_hdag = [res for res in results if res[1] > 0]
                num_nodes = len(results_in_hdag)

            with open(real_data_results_path, "rb") as f:
                results = pickle.load(f)
                results_in_hdag = [res for res in results if res[1] > 0]
                real_data_num_nodes = len(results_in_hdag)

            with open(sim_info_path, "r") as f:
                tree_stats = json.load(f)
                toi_score = int(tree_stats['max_score_top'])
                toi_num_nodes = int(tree_stats['num_nodes'])
                toi_num_leaves = int(tree_stats['num_leaves'])
                toi_num_internal = toi_num_nodes - toi_num_leaves

            with open(support_log_path, "r") as f:
                line = f.read().split("=>")[-1].strip()
                assert "DAG contains " in line
                num_trees = int(line.split(" ")[2])

            
            trial_list.append((clade_name, prob * rate, trial, num_trees, num_nodes, toi_num_internal, \
                toi_num_leaves, real_num_leaves, toi_num_leaves / real_num_leaves, toi_score, mp_score, real_mp_score, \
                (num_nodes / toi_num_leaves) / (real_data_num_nodes / real_num_leaves), \
                toi_score / mp_score, real_data_num_nodes, mp_score / real_mp_score))

        except Exception as e:
            print(e)
            continue

    # trial_list.sort(key=lambda el: el[3], reverse=True)
    # trial_list.sort(key=lambda el: el[4], reverse=True)

    df = pd.DataFrame(trial_list, columns =["clade", "sim_type", "trial", \
            "num_trees", "num_nodes", "toi_num_internal", "toi_num_leaves", "real_num_leaves", \
            "prop_new_leaves", "toi_score", "mp_score", "real_mp_score", \
            "prop_new_nodes", "suboptimality", "real_data_num_nodes", "prop_mp"])
    print(df.to_string())

    fig, axes = plt.subplots()

    x_col = "sim_type"
    sns.boxplot(x=x_col, y='prop_new_nodes', data=df, ax=axes)

    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Difference in Parsimony Diversity (sim / real data)')
    axes.set_yscale("log", base=2)
    # plt.xticks(rotation=15)
    # axes.set_xscale("log")

    plt.savefig(data_path + f"/figures/real_comparison_pars_div_by_{x_col}_boxplot_{len(clade_names)}_clades.png")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x=x_col, y='num_nodes', data=df, ax=axes)

    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Parsimony Diversity (Number of DAG nodes)')
    # plt.xticks(rotation=15)

    plt.savefig(data_path + f"/figures/pars_div_by_{x_col}_boxplot_{len(clade_names)}_clades.png")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x=x_col, y='suboptimality', data=df, ax=axes)

    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Suboptimality (TOI score / MP score)')
    # plt.xticks(rotation=15)

    plt.savefig(data_path + f"/figures/suboptimality_by_{x_col}_boxplot_{len(clade_names)}_clades.png")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x=x_col, y='prop_mp', data=df, ax=axes)

    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Simulated MP score / Real MP score')
    # plt.xticks(rotation=15)

    plt.savefig(data_path + f"/figures/prop_mp_{x_col}_boxplot_{len(clade_names)}_clades.png")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x=x_col, y='prop_new_leaves', data=df, ax=axes)

    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Simulated Numn Leaves / Real Num Leaves')
    # plt.xticks(rotation=15)

    plt.savefig(data_path + f"/figures/prop_new_leaves_{x_col}_boxplot_{len(clade_names)}_clades.png")
    plt.clf()


@click.command("hmut_coverage")
@click.option("-d", "--data-path", help="Path to data directory.")
@click.option("-o", "--out-dir", help='output path to store figures/tables in.')
@click.option("-m", "--method", default='historydag')
@click.option('--support-removal-threshold', '-t', default=0.01)
def hmut_coverage(data_path, out_dir, method, support_removal_threshold):
    """
    E.g., 
python support_pipeline/plotting.py hmut_coverage \
-d "/fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation" \
-o "/fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation/figures" \
-m mrbayes
    """
    clade_names = [
        "AY.34.2",
        "AY.108",
        "AZ.3",
        "P.1.7"
    ]

    hmut_params = [
        (0.01, 200),
        (0.001, 200),
        (0.0001, 200),
        (0.01, 20),
        (0.001, 20),
        (0.0001, 20),
        # (0.01, 100),
        # (0.001, 100),
        # (0.0001, 100),
        # (0.01, 10),
        # (0.001, 10),
        # (0.0001, 10)
    ]

    # hmut_params = [h for h in hmut_params if h[0] * h[1] >= 0.1]

    trials = [1, 2, 3]

    param_trial_list = [(el1, el2, el3) for el1 in clade_names for el2 in hmut_params for el3 in trials]

    hmut2results = {}

    # Iterate over every combination of hmut param and trial
    for i, (clade_name, (prob, rate), trial) in enumerate(param_trial_list):
        base_path = f"{data_path}/{clade_name}/{prob}_{rate}_JC/{trial}"
        results_path = base_path + f"/results/{method}/results.pkl"

        if i % 2 == 0:
            print(f"{i} / {len(param_trial_list)}")

        try:
            with open(results_path, "rb") as f:
                results = pickle.load(f)
                results = [result for result in results if len(result[0]) > 1 and result[1] > support_removal_threshold]
        except:
            print(f"No results at {base_path}")

        if prob * rate not in hmut2results:
            hmut2results[prob * rate] = []
        hmut2results[prob * rate].extend(results)

    print("Creating CA plot...")
    window_size=100

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[0.75, 0.25])
    fig.set_size_inches(7, 9)
    ax1.set_title(f"Coverage analysis for different hypermutation settings")
    ax1.set_ylabel(f"Empirical Probability")
    ax2.set_yscale("log")
    ax2.set_xlabel("Estimated Support")

    with open(f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation/figures/{method}_coverage.pkl", "wb") as f:
        pickle.dump(hmut2results, f)
    with open(f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/hypermutation/figures/{method}_coverage.pkl", "rb") as f:
        hmut2results = pickle.load(f)

    hist_list = []
    for hmut, results_full in hmut2results.items():
        print(hmut)
        if method == 'historydag':
            window_size = int(len(results_full) / 20)
        x, y, _, _ = sliding_window_plot(results_full, std_dev=True, window_size=window_size)
        ax1.plot(x, y, label=f"{hmut}")
        hist_list.append(x)
    
    ax2.hist(hist_list)
    ax1.plot([0, 1], [0, 1], color="blue", label="Ideal", linestyle="dashed")
    ax1.legend()
    
    fig.savefig(out_dir + f"/{method}_hmut_coverage.png")
    fig.clf()


#################################################
### Simulation Model Experiments
#################################################

@click.command("plot_pars_distribution")
@click.option("-d", "--data-path", help="Path to data directory.")
@click.option("-i", "--input-name", default="mrbayes/mrbayes-pars-distr.pkl", \
              help="Path to sample file to concatenate after clade, trial, sim_type.")
def plot_pars_distribution(data_path, input_name):
    """
    Plot the distribution of parsimony scores for posterior given by samples at data_path
    """
    clade_names = [
        "AY.34.2",
        "AY.108",
        "AZ.3",
        "P.1.7",
        "AY.74",
        "A.2.5",
        "AY.87",
        "B.1.1.10"
    ]

    sim_params = [
        # "jc",
        # "jc_rv",
        # "jc_rv_inv",
        # "jc_rv_inv_hmut",
        # "gtr",
        # "unrest_hmut"
        "gamma_10_hmut_50"
    ]

    trials = range(1, 26)

    param_trial_list = [(el1, el2, el3) for el1 in clade_names for el2 in sim_params for el3 in trials]

    clade_sim_2_change = {}

    # Iterate over every combination of hmut param and trial
    for clade_name, sim_type, trial in param_trial_list:
        if trial == 1:
            print(clade_name, sim_type)
        base_path = f"{data_path}/{clade_name}/{sim_type}/{trial}"

        pars_distr_path = base_path + f"/results/{input_name}"
        sim_info_path = base_path + f"/simulation/tree_stats.json"

        fig_path = base_path + f"/figures/diff_historydag_mut_rates"    # TODO: May need to adjust this

        try:
            with open(pars_distr_path, "rb") as f:
                pars_distr = pickle.load(f)

            with open(sim_info_path, "r") as f:
                tree_stats = json.load(f)
                toi_score = int(tree_stats['max_score_top'])

        except Exception as e:
            print(e)
            continue
        
        # TODO:
        # - plot distribution
        # - aggregate

        if not os.path.isdir(fig_path):
            os.makedirs(fig_path)

        total_count = sum(pars_distr.values())
        width = 1
        indexes = np.array(range(min(int(min(pars_distr.keys())), toi_score), max(int(max(pars_distr.keys())), toi_score)+1))
        values = []
        labels = []
        relative_change_in_pars = []
        for i in indexes:
            labels.append(i)
            if i in pars_distr:
                values.append(pars_distr[i] / total_count)
            else:
                values.append(0)

            relative_change_in_pars.append(((i - toi_score)/toi_score, values[-1] * total_count))
        
        if (clade_name, sim_type) not in clade_sim_2_change:
            clade_sim_2_change[(clade_name, sim_type)] = []
        clade_sim_2_change[(clade_name, sim_type)].extend(relative_change_in_pars)

        f = plt.figure(figsize=(6.4, 5.8))
        plt.xticks(indexes, labels, rotation=45)
        plt.bar(indexes, values, width)

        plt.axvline(toi_score, color='red', linestyle='dashed', linewidth=1)

        plt.xlabel("Parsimony Score")
        plt.ylabel("Posterior")
        plt.title(f"Parsimony Posterior for {clade_name} {sim_type}")
        plt.savefig(f"{fig_path}/parsimony_distribution.png")
        plt.clf()
        plt.close()

    for clade_name, sim_type, trial in param_trial_list:
        if trial != 1:
            continue

        base_path = f"{data_path}/{clade_name}/{sim_type}"
        fig_path = base_path + f"/figures/diff_historydag_mut_rates"
        if not os.path.isdir(fig_path):
            os.makedirs(fig_path)

        relative_change_list = clade_sim_2_change[(clade_name, sim_type)]
        # print(relative_change_list)
        
        # Combine relative changes if they're the same value
        # for change in relative_change_list:
        relative_change = {}
        for change, count in relative_change_list:
            if change not in relative_change:
                relative_change[change] = 0
            relative_change[change] += count/(len(trials) * total_count)
        
        num_bins=10
        min_change = min(relative_change.keys()) * 0.999
        max_change = max(relative_change.keys()) * 1.001
        spaces = np.linspace(min_change, max_change, num=num_bins+1)
        bar = {}
        for i, start in enumerate(spaces[:-1]):
            val = 0
            for change in relative_change.keys():
                if change >= spaces[i] and change < spaces[i+1]:
                    val += relative_change[change]
            bar[start] = val

        # print(bar)

        f = plt.figure(figsize=(6.4, 7))
        plt.bar(bar.keys(), bar.values(), width=spaces[1]-spaces[0], align="edge")
        plt.xticks(rotation=30)

        plt.axvline(0, color='red', linestyle='dashed', linewidth=2)

        plt.xlabel("Relative Change in Parsimony Score (Sample - Sim) / Sim")
        plt.ylabel("Posterior")
        plt.title(f"Posterior on Relative Change in Parsimony {clade_name} {sim_type}")
        plt.savefig(f"{fig_path}/parsimony_distribution_relative_change.png")
        plt.clf()
        plt.close()
    
    # NOTE: For plotting multiple simulation types
    # sim2idx = {sim: (i//3, i%3) for i, sim in enumerate(sim_params)}
    # for curr_clade in clade_names:
    #     print(curr_clade)
    #     fig, axes = plt.subplots(2, 3)
    #     fig.set_size_inches(8, 6)
    #     for clade_name, sim_type, trial in param_trial_list:
    #         if trial != 1 or curr_clade != clade_name:
    #             continue
    #         print("\t", sim_type)
    #         base_path = f"{data_path}/{clade_name}"
    #         fig_path = base_path + f"/figures"
    #         if not os.path.isdir(fig_path):
    #             os.makedirs(fig_path)

    #         relative_change_list = clade_sim_2_change[(clade_name, sim_type)]
    #         print(len(relative_change_list))
            
    #         # Combine relative changes if they're the same value
    #         # for change in relative_change_list:
    #         relative_change = {}
    #         for change, count in relative_change_list:
    #             if change not in relative_change:
    #                 relative_change[change] = 0
    #             relative_change[change] += count/(5 * total_count)
            
    #         num_bins=10
    #         min_change = min(relative_change.keys()) * 0.999
    #         max_change = max(relative_change.keys()) * 1.001
    #         spaces = np.linspace(min_change, max_change, num=num_bins+1)
    #         bar = {}
    #         for i, start in enumerate(spaces[:-1]):
    #             val = 0
    #             for change in relative_change.keys():
    #                 if change >= spaces[i] and change < spaces[i+1]:
    #                     val += relative_change[change]
    #             bar[start] = val

    #         # plt.xticks(rotation=30)
    #         ax = axes[sim2idx[sim_type]]
    #         ax.bar(bar.keys(), bar.values(), width=spaces[1]-spaces[0], align="edge")
    #         ax.axvline(0, color='red', linestyle='dashed', linewidth=2)
    #         ax.set_title(sim_type)

    #     fig.supxlabel("Relative Change in Parsimony Score (Sample - Sim) / Sim")
    #     fig.supylabel("MrBayes Posterior")
    #     fig.suptitle(f"MB Posterior on Relative Change in Parsimony {curr_clade}")
    #     fig.tight_layout()
    #     plt.savefig(f"{fig_path}/parsimony_distribution_relative_change.png")
    #     plt.clf()
    #     plt.close()



@click.command("plot_sim_model_stats")
@click.option("-d", "--data-path", help="Path to data directory.")
def plot_sim_model_stats(data_path):
    """
    Gathers statistics about hypermutation simulations and creates figures summarizing the results.
    """

    clade_names = [
        "AY.34.2",
        "AY.108",
        "AZ.3",
        "P.1.7"
    ]

    hmut_params = [
        # "jc",
        # "jc_rv",
        # "jc_rv_inv",
        # "jc_rv_inv_hmut",
        # "gtr",
        "unrest_hmut"
    ]

    trials = [1, 2, 3, 4, 5]

    param_trial_list = [(el1, el2, el3) for el1 in clade_names for el2 in hmut_params for el3 in trials]

    trial_list = []

    # Iterate over every combination of hmut param and trial
    for clade_name, sim_type, trial in param_trial_list:
        base_path = f"{data_path}/{clade_name}/{sim_type}/{trial}"

        log_path = base_path + "/results/historydag/opt_info/optimization_log_complete/logfile.csv"
        sim_info_path = base_path + f"/simulation/tree_stats.json"
        results_path = base_path + f"/results/historydag/results.pkl"
        support_log_path = base_path + f"/results/historydag/inference_historydag.log"
        real_data_results_path = data_path + f"/../{clade_name}_/results/historydag/results.pkl"

        try:
            with open(log_path, "r") as f:
                mp_score = int(f.readlines()[-1].split("\t")[4])

            with open(results_path, "rb") as f:
                results = pickle.load(f)
                results_in_hdag = [res for res in results if res[1] > 0]
                num_nodes = len(results_in_hdag)

            with open(real_data_results_path, "rb") as f:
                results = pickle.load(f)
                results_in_hdag = [res for res in results if res[1] > 0]
                real_data_num_nodes = len(results_in_hdag)

            with open(sim_info_path, "r") as f:
                tree_stats = json.load(f)
                toi_score = int(tree_stats['max_score_top'])
                toi_num_nodes = int(tree_stats['num_nodes'])
                toi_num_leaves = int(tree_stats['num_leaves'])
                toi_num_internal = toi_num_nodes - toi_num_leaves

            with open(support_log_path, "r") as f:
                line = f.read().split("=>")[-1].strip()
                assert "DAG contains " in line
                num_trees = int(line.split(" ")[2])

            if toi_score < mp_score:
                print(f"{toi_score} and {mp_score}: Did not find MP trees in", clade_name, sim_type, trial)
                continue
            
            trial_list.append((clade_name, sim_type, trial, \
                num_trees, num_nodes, toi_num_internal, toi_score, mp_score, \
                num_nodes / real_data_num_nodes, toi_score / mp_score, real_data_num_nodes))

        except Exception as e:
            print(e)
            continue

    trial_list.sort(key=lambda el: el[3], reverse=True)
    trial_list.sort(key=lambda el: el[4], reverse=True)

    df = pd.DataFrame(trial_list, columns =["clade", "sim_type", "trial", \
            "num_trees", "num_nodes", "toi_num_internal", "toi_score", "mp_score", \
            "prop_new_nodes", "suboptimality", "real_data_num_nodes"])
    print(df)

    fig, axes = plt.subplots()
    sns.boxplot(x='sim_type',y='prop_new_nodes', data=df, ax=axes)

    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Difference in Parsimony Diversity (sim / real data)')
    axes.set_yscale("log", base=2)
    plt.xticks(rotation=15)
    # axes.set_xscale("log")

    plt.savefig(data_path + "/figures/real_comparison_pars_div_boxplot.png")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x='sim_type',y='num_nodes', data=df, ax=axes)

    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Parsimony Diversity (Number of DAG nodes)')
    plt.xticks(rotation=15)

    plt.savefig(data_path + "/figures/pars_div_boxplot.png")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x='sim_type',y='suboptimality', data=df, ax=axes)

    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Suboptimality (TOI score / MP score)')
    plt.xticks(rotation=15)

    plt.savefig(data_path + "/figures/suboptimality_boxplot.png")
    plt.clf()

@click.command("sim_model_coverage")
@click.option("-d", "--data-path", help="Path to data directory.")
@click.option("-o", "--out-dir", help='output path to store figures/tables in.')
@click.option("-m", "--method", default='historydag')
@click.option('--support-removal-threshold', '-t', default=0.01)
def sim_model_coverage(data_path, out_dir, method, support_removal_threshold):
    """
    E.g., 
python support_pipeline/plotting.py sim_model_coverage \
-d "/fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models" \
-o "/fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models/figures" \
-m historydag
    """
    # clade_names = [
    #     "AY.34.2",
    #     "AY.108",
    #     "AZ.3",
    #     "P.1.7"
    # ]
    # TODO: Refactor hmut --> sim_models
    # hmut_params = [
    #     "jc",
    #     "jc_rv",
    #     "jc_rv_inv",
    #     "jc_rv_inv_hmut",
    #     "gtr",
    #     "unrest_hmut"
    # ]

    clade_names = [
        "AY.34.2",
        "AY.108",
        "AZ.3",
        "P.1.7",
        "A.2.5",
        "AY.74",
        "AY.87",
        "B.1.1.10"
    ]
    hmut_params = [
        "unrest_hmut_150",
        "unrest_hmut_170",
        "unrest_hmut_200",
        "unrest_hmut_220"
    ]

    trials = [1, 2, 3, 4, 5]

    param_trial_list = [(el1, el2, el3) for el1 in clade_names for el2 in hmut_params for el3 in trials]

    hmut2results = {}

    # TODO: Uncomment
    # Iterate over every combination of hmut param and trial
    for i, (clade_name, sim_type, trial) in enumerate(param_trial_list):
        base_path = f"{data_path}/{clade_name}/{sim_type}/{trial}"
        results_path = base_path + f"/results/{method}/results.pkl"

        if i % 2 == 0:
            print(f"{i} / {len(param_trial_list)}")

        try:
            with open(results_path, "rb") as f:
                results = pickle.load(f)
                results = [result for result in results if len(result[0]) > 1 and result[1] > support_removal_threshold]
        except:
            print(f"No results at {base_path}")

        if sim_type not in hmut2results:
            hmut2results[sim_type] = []
        hmut2results[sim_type].extend(results)
    with open(f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models/figures/{method}_coverage.pkl", "wb") as f:
        pickle.dump(hmut2results, f)
    with open(f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/sim_models/figures/{method}_coverage.pkl", "rb") as f:
        hmut2results = pickle.load(f)

    print("Creating CA plot...")
    window_size=100

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[0.75, 0.25])
    fig.set_size_inches(7, 9)
    ax1.set_title(f"Coverage analysis for different simulation settings")
    ax1.set_ylabel(f"Empirical Probability")
    ax2.set_yscale("log")
    ax2.set_xlabel("Estimated Support")


    hist_list = []
    for hmut, results_full in hmut2results.items():
        print(hmut)
        if True: # method == 'historydag':
            window_size = int(len(results_full) / 20)
        x, y, _, _ = sliding_window_plot(results_full, std_dev=True, window_size=window_size)
        ax1.plot(x, y, label=f"{hmut}")
        hist_list.append(x)
    
    ax2.hist(hist_list)
    ax1.plot([0, 1], [0, 1], color="blue", label="Ideal", linestyle="dashed")
    ax1.legend()
    
    fig.savefig(out_dir + f"/{method}_sim_models_coverage.png")
    fig.clf()


@click.command("plot_sim_stats")
@click.option("-d", "--data-path", help="Path to data directory.")
def plot_sim_stats(data_path):
    """
    Gather info about complex simulations with a variety of hypermutation scores and plot figures
    summarizing that information.
    """

    num_trials = 5
    clade_names = [
        "AY.34.2",
        "AY.108",
        "AZ.3",
        "P.1.7",
        "A.2.5",
        "AY.74",
        "AY.87",
        "B.1.1.10"
    ]

    hmut_params = [
        # hmut varies, gamma = 0
        "unrest_hmut_0",
        "unrest_hmut_50_10",
        "unrest_hmut_100_10",
        "unrest_hmut_150_10",
        # hmut = 0, gamma varies
        "gamma_rv_10",
        "gamma_rv_20",
        # gamma and hmut vary
        "gamma_5_hmut_50",
        "gamma_5_hmut_100",
        "gamma_5_hmut_150",
        "gamma_5_hmut_200",
        "gamma_10_hmut_50",
        "gamma_10_hmut_100",
        "gamma_10_hmut_150",
        "gamma_10_hmut_200",
        "gamma_20_hmut_50",
        "gamma_20_hmut_100",
        "gamma_20_hmut_150",
        "gamma_20_hmut_200",
        "gamma_30_hmut_50",
        "gamma_30_hmut_100",
        "gamma_30_hmut_150",
        "gamma_30_hmut_200",
    ]

    trials = range(1, num_trials+1)

    param_trial_list = [(el1, el2, el3) for el1 in clade_names for el2 in hmut_params for el3 in trials]

    trial_list = []

    # Iterate over every combination of hmut param and trial
    for clade_name, sim_type, trial in param_trial_list:
        base_path = f"{data_path}/{clade_name}/{sim_type}/{trial}"

        log_path = base_path + "/results/historydag/opt_info/optimization_log_complete/logfile.csv"
        sim_info_path = base_path + f"/simulation/tree_stats.json"
        results_path = base_path + f"/results/historydag/results.pkl"
        support_log_path = base_path + f"/results/historydag/inference_historydag.log"
        real_data_results_path = data_path + f"/../real_data/{clade_name}/results/historydag/results.pkl"
        real_data_hdag_log_path = data_path + f"/../real_data/{clade_name}/results/historydag/opt_info/optimization_log_complete/logfile.csv"
        real_data_fasta_path = data_path + f"/../real_data/{clade_name}/reconstructed_seqs.fasta"

        try:
            with open(real_data_hdag_log_path, "r") as f:
                real_mp_score = int(f.readlines()[-1].split("\t")[4])

            with open(real_data_fasta_path, "r") as f:
                real_num_leaves = int(len(f.readlines())/2)

            with open(log_path, "r") as f:
                mp_score = int(f.readlines()[-1].split("\t")[4])

            with open(results_path, "rb") as f:
                results = pickle.load(f)
                results_in_hdag = [res for res in results if res[1] > 0]
                num_nodes = len(results_in_hdag)

            with open(real_data_results_path, "rb") as f:
                results = pickle.load(f)
                results_in_hdag = [res for res in results if res[1] > 0]
                real_data_num_nodes = len(results_in_hdag)

            with open(sim_info_path, "r") as f:
                tree_stats = json.load(f)
                toi_score = int(tree_stats['pars_score'])
                max_score_for_top = int(tree_stats['max_score_top'])
                toi_num_nodes = int(tree_stats['num_nodes'])
                toi_num_leaves = int(tree_stats['num_leaves'])
                toi_num_internal = toi_num_nodes - toi_num_leaves

            with open(support_log_path, "r") as f:
                line = f.read().split("=>")[-1].strip()
                assert "DAG contains " in line
                num_trees = int(line.split(" ")[2])

            if toi_score < mp_score:
                print(f"{toi_score} and {mp_score}: Did not find MP trees in", clade_name, sim_type, trial)
                continue
            
            # if sim_type.split("_")[2] != '0' and sim_type.split("_")[3] == '5':
            #     continue

            # NOTE: May have to change these depending on how the simulations are named
            # if sim_type == "unrest":
            #     sim_type = "No RV"
            # elif sim_type[0:3] == "gam":
            #     arr = sim_type.split("_")
            #     sim_type = arr[0] + " " + str(int(arr[2]) / 100)
            # else:
            #     arr = sim_type.split("_")
            #     sim_type = arr[2] + " " + arr[3]

            gamma = None
            hmut = None
            if sim_type.startswith("unrest"):
                gamma = 0
                hmut = int(sim_type.split("_")[2])
            elif sim_type.startswith("gamma_rv"):
                gamma = int(sim_type.split("_")[2]) / 100
                hmut = 0
            else:
                gamma = int(sim_type.split("_")[1]) / 100
                hmut = int(sim_type.split("_")[3])
            
            sim_type = (hmut, gamma)

            trial_list.append((clade_name, sim_type, trial, num_trees, num_nodes, toi_num_internal, \
                toi_num_leaves, real_num_leaves, toi_num_leaves / real_num_leaves, toi_score, mp_score, real_mp_score, \
                (num_nodes / toi_num_leaves) / (real_data_num_nodes / real_num_leaves), \
                toi_score / mp_score, real_data_num_nodes, (mp_score / toi_num_leaves) / (real_mp_score/ real_num_leaves), \
                ((num_nodes / toi_num_leaves) - (real_data_num_nodes / real_num_leaves)) / (real_data_num_nodes / real_num_leaves),
                max_score_for_top / toi_score))

        except Exception as e:
            print(e)
            continue
    
    trial_list.sort(reverse=True, key=lambda el: el[12])
    df = pd.DataFrame(trial_list, columns =["clade", "sim_type", "trial", \
                "num_trees", "num_nodes", "toi_num_internal", "toi_num_leaves", "real_num_leaves", \
                "prop_new_leaves", "toi_score", "mp_score", "real_mp_score", "prop_new_nodes", \
                "suboptimality", "real_data_num_nodes", "prop_mp", "rel_change_new_nodes", \
                "consistency_index"])
    # print(df.to_string())

    sns.set_theme()
    sns.set(font_scale=1.2)
    sns.set_style("white")
    x_col = "sim_type"

    fig, axes = plt.subplots()
    fig.set_size_inches(10, 10)
    sns.boxplot(x=x_col, y='prop_new_nodes', data=df, ax=axes)
    axes.yaxis.grid(True)
    axes.set_xlabel('Parameter Values')
    axes.set_ylabel('CIPD')
    axes.set_yscale("log", base=2)
    plt.xticks(rotation=45)
    plt.yticks([2 ** i for i in range(-2, 8, 2)])
    # plt.title("Difference in Parsimony Diversity")
    plt.savefig(data_path + f"/figures/real_comparison_pars_div_by_{x_col}_boxplot_{len(clade_names)}_clades.pdf")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x=x_col, y='rel_change_new_nodes', data=df, ax=axes)
    axes.yaxis.grid(True)
    axes.set_xlabel('Hypermutation Rate')
    axes.set_ylabel('Relative Change in Parsimony Diversity')
    axes.set_ylim(top = 7)
    # axes.set_yscale("log", base=2)
    plt.savefig(data_path + f"/figures/rel_change_pars_div_{x_col}_boxplot_{len(clade_names)}_clades.png")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x=x_col, y='num_nodes', data=df, ax=axes)
    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Parsimony Diversity (Number of DAG nodes)')
    plt.xticks(rotation=15)
    plt.savefig(data_path + f"/figures/pars_div_by_{x_col}_boxplot_{len(clade_names)}_clades.png")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x=x_col, y='suboptimality', data=df, ax=axes)
    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Suboptimality (TOI score / MP score)')
    plt.xticks(rotation=15)
    plt.savefig(data_path + f"/figures/suboptimality_by_{x_col}_boxplot_{len(clade_names)}_clades.png")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x=x_col, y='prop_mp', data=df, ax=axes)
    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Simulated MP score / Real MP score')
    plt.xticks(rotation=15)
    plt.savefig(data_path + f"/figures/prop_mp_{x_col}_boxplot_{len(clade_names)}_clades.png")
    plt.clf()

    fig, axes = plt.subplots()
    sns.boxplot(x=x_col, y='prop_new_leaves', data=df, ax=axes)
    axes.yaxis.grid(True)
    axes.set_xlabel('Simulation Type')
    axes.set_ylabel('Simulated Numn Leaves / Real Num Leaves')
    # plt.xticks(rotation=15)

    plt.savefig(data_path + f"/figures/prop_new_leaves_{x_col}_boxplot_{len(clade_names)}_clades.png")
    plt.clf()

    heat_map_avg = np.zeros((4, 4))
    heat_map_med = np.zeros((4, 4))
    xs = df['sim_type'].to_numpy()
    ys = df['prop_new_nodes'].to_numpy()
    alpha_vals = [0, 0.05, 0.1, 0.2, 0.3] # [0, 0.3, 0.2, 0.1]
    # alpha_vals.sort(reverse=True)
    hmut_vals = [0, 50, 100, 150, 200]
    # alpha_vals.sort(reverse=True)
    for a in alpha_vals:
        for h in hmut_vals:
            key = (h, a)
            trials = []
            for x, y in zip(xs, ys):
                if x == key:
                    trials.append(y)
            if len(trials) != 0:
                avg = sum(trials) / len(trials)
                trial_arr = np.array(trials)
                med = np.median(trial_arr)
            else:
                med = np.nan
                avg = np.nan
            heat_map_med[alpha_vals.index(a), hmut_vals.index(h)] = med
            heat_map_avg[alpha_vals.index(a), hmut_vals.index(h)] = avg
    
    ax = sns.heatmap(heat_map_avg, linewidth=0.5, yticklabels=alpha_vals, xticklabels=hmut_vals,
                     annot=True, vmin=0.5)
    ax.set_title("Average Parsimony Diversity")
    ax.set_xlabel("Hypermutation Rate")
    ax.set_ylabel("Alpha")
    plt.savefig(data_path + f"/figures/PD_heat_map_avg_{len(clade_names)}_clades.png")
    plt.clf()

    ax = sns.heatmap(heat_map_med, linewidth=0.5, yticklabels=alpha_vals, xticklabels=hmut_vals,
                     annot=True, vmin=0.5, cbar_kws={'label': 'Corrected PD Increase'})
    # ax.set_title("Median Parsimony Diversity")
    ax.set_xlabel("Hypermutation Rate")
    ax.set_ylabel("Alpha")
    plt.savefig(data_path + f"/figures/PD_heat_map_med_{len(clade_names)}_clades.pdf")
    plt.clf()

    # TODO: Put these in different methods so you can comment them out easily
    # y_col = 'prop_new_nodes'
    # xs = df['sim_type'].to_numpy()
    # ys = df[y_col].to_numpy()
    # clades = df['clade'].to_numpy()

    # import matplotlib.cm as cm
    # clade_colors = cm.tab10(np.linspace(0, 1, len(clade_names)))
    # clade2color = {clade: color for clade, color in zip(clade_names, clade_colors)}
    # # colors = [clade2color[clade] for clade in clades]
    # # plt.scatter(xs, ys, color=colors, alpha=0.4)

    # clade2idx = {clade: np.where(clades == clade)[0] for clade in clade_names}

    # plots = []
    # for clade in clade_names:
    #     x_vals = [150, 200] #[0, 50, 100, 150, 200, 250]
    #     y_vals = []
    #     for hmut in x_vals:
    #         idxs = np.where((clades == clade) & (xs == hmut))[0]
    #         y_vals.append(np.mean(ys[idxs]))
    #     plots.append(plt.scatter(x_vals, y_vals, color=clade2color[clade], alpha=0.7))

    # plt.legend(plots,
    #        clade_names,
    #        scatterpoints=1,
    #     #    loc='lower left',
    #        ncol=2)

    # plt.plot([min(x_vals) * 0.99, max(x_vals) * 1.01], [1, 1], color="red", linestyle='dashed')

    # plt.xlabel('Hypermutation Rate')
    # plt.ylabel('Corrected Proportion New Nodes (sim / real)')
    # # plt.yscale('log', base=2)
    # plt.savefig(data_path + f"/figures/{y_col}_scatter_{len(clade_names)}_clades.png")
    # plt.clf()

    # Parsimony Diversity vs Consistency Index scatter
    # fig = plt.figure()
    # fig.set_size_inches(7.5, 5)
    # xs = df['prop_new_nodes'].to_numpy()
    # ys = df['consistency_index'].to_numpy()
    # clades = df['clade'].to_numpy()

    # plots = []
    # for clade in clade_names:
    #     idxs = np.where((clades == clade))[0]
    #     x_vals = xs[idxs]
    #     y_vals = ys[idxs]
    #     plots.append(plt.scatter(x_vals, y_vals, color=clade2color[clade], alpha=0.7))

    # plt.legend(plots,
    #        clade_names,
    #        scatterpoints=1,
    #        ncol=2)

    # plt.xlabel('Corrected Parsimony Diversity')
    # plt.ylabel('Consistency Index (Minimum # muts / actual # muts)')
    # plt.xscale('log', base=2)
    # plt.savefig(data_path + f"/figures/PD_CI_scatter_{len(clade_names)}_clades.png")
    # plt.clf()


#################################################
### Suboptimal Structure Experiments
#################################################
    
@click.command("plot_sub_stats")
@click.option("-d", "--data-path", help="Path to data directory.")
def plot_sub_stats(data_path):
    """
    Gather info about complex simulations with a variety of hypermutation scores and plot figures
    summarizing that information.
    """

    num_trials = 25
    clade_names = [
        "AY.34.2",
        "AY.108",
        "AZ.3",
        "P.1.7",
        "A.2.5",
        "AY.74",
        "AY.87",
        "B.1.1.10"
    ]
    hmut_params = [
        "gamma_10_hmut_50"
    ]
    trials = range(1, num_trials+1)
    param_trial_list = [(el1, el2, el3) for el1 in clade_names for el2 in hmut_params for el3 in trials]

    def parse_dups(f):
        line = f.readline()
        while not line.startswith("===") and line:
            line = f.readline()
            # print(line)
        # print(f.readline().strip().split("\t"))
        diffs = int(f.readline().strip().split("\t")[1])
        dups = int(f.readline().strip().split("\t")[1])
        return diffs, dups

    diff_counts, dup_counts, name_trial = [], [], []

    # Iterate over every combination of hmut param and trial
    for clade_name, sim_type, trial in param_trial_list:
        base_path = f"{data_path}/{clade_name}/{trial}"
        log_path = base_path + "/compare_trial.log"

        try:
            with open(log_path, "r") as f:
                diffs, dups = parse_dups(f)
        except Exception as e:

            print(f"MISSING: {clade_name} {trial} {e}")
            continue

        diff_counts.append(diffs)
        dup_counts.append(dups)
        name_trial.append((clade_name, trial))

    print(f"{sum(dup_counts)} out of {sum(diff_counts)} ({sum(dup_counts) / sum(diff_counts) * 100:.2f}%) are due to PCM")

    proportion_dup = np.array([dup / diff if diff > 0 else 1 for dup, diff in zip(dup_counts, diff_counts)])
    srtd_idxs = np.argsort(proportion_dup)
    print("Top 10 worst trials are")
    for i in srtd_idxs[:10]:
        print("\t", name_trial[i][0], name_trial[i][1], f"\t{dup_counts[i]} / {diff_counts[i]} \t= {proportion_dup[i]}")

    print()
    print("Proportion diffs due to PCM per clade")
    for curr_clade in clade_names:
        idxs = [i for i, (clade, _) in enumerate(name_trial) if clade == curr_clade]
        total_diffs = 0
        total_dups = 0
        for i in idxs:
            total_diffs += diff_counts[i]
            total_dups += dup_counts[i]

        print(f"\t{curr_clade}\t{total_dups} / {total_diffs} \t = {total_dups / total_diffs}")

    clades = np.array([clade for clade, _ in name_trial])
    # Sort
    proportion_dup = proportion_dup[srtd_idxs]
    clades = clades[srtd_idxs]
    
    sns.set_theme(style="whitegrid", font_scale=1.25)
    df = pd.DataFrame({"prop_dup": proportion_dup, "clade": clades})
    fig, axes = plt.subplots()
    
    # fig.set_size_inches((8.5, 4.8))

    
    sns.swarmplot(x="clade", y='prop_dup', hue='clade', data=df, ax=axes, legend=False, size=4)
    df_means = df.groupby("clade")["prop_dup"].agg("mean").reset_index()
    xlim = axes.get_xlim()
    ylim = axes.get_ylim()
    sns.scatterplot(x="clade", y="prop_dup", marker='X', color='black', s=130, zorder=3, ax=axes, legend=False, data=df_means)
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    plt.xticks(rotation=15)

    plt.ylabel("Proportion Duplicate Nodes")
    plt.xlabel("UShER Clade")
    plt.savefig(f"{data_path}/figures/proportion_duplicate_swarm_per_clade.pdf")
    plt.clf()

    plt.hist(proportion_dup)
    plt.savefig(f"{data_path}/figures/proportion_duplicate_histogram.pdf")
    plt.clf()

#################################################
### Diffused DAG Sampling Experiments
#################################################

@click.command("compare_mutation_probs")
@click.option("-d", "--data-path", help="Path to data directory.")
@click.option("-i", "--input-name",  help="Path to sample file to concatenate after clade, trial, sim_type.")
def compare_mutation_probs(data_path, input_name):
    """
    Plot the distribution of parsimony scores for posterior given by samples at data_path
    """
    clade_names = [
        "AY.34.2",
        "AY.108",
        "AZ.3",
        "P.1.7",
        "AY.74",
        "A.2.5",
        "AY.87",
        "B.1.1.10"
    ]
    sim_params = [
        "gamma_10_hmut_50"
    ]
    trials = range(1, 26)
    param_trial_list = [(el1, el2, el3) for el1 in clade_names for el2 in sim_params for el3 in trials]

    # Iterate over every combination of hmut param and trial
    for clade_name, sim_type, trial in param_trial_list:
        if trial == 1:
            print(clade_name, sim_type)
        base_path = f"{data_path}/{clade_name}/{sim_type}/{trial}"

        probs_dict_path = base_path + f"/results/{input_name}"

        fig_path = base_path + f"/figures/diff_historydag"

        try:
            with open(probs_dict_path, "rb") as f:
                probs_dict = pickle.load(f)

        except Exception as e:
            print(e)
            continue

        if not os.path.isdir(fig_path):
            os.makedirs(fig_path)


        f = plt.figure(figsize=(6.4, 5.8))
        bins = np.geomspace(1e-5, 1, 10)
        # bins = np.linspace(0, 1, 20)
        hist_list = []
        label_list = []
        for name, probs in probs_dict.items():
            hist_list.append(probs)
            label_list.append(name)
        plt.hist(hist_list, bins, label=label_list)

        plt.legend()
        plt.xlabel("Probability of Mutation")
        plt.ylabel("Counts")
        # plt.xscale("log")
        plt.yscale("log")
        plt.title(f"Comparsion of Different Mutation Probability Functions {clade_name} {sim_type}")
        plt.savefig(f"{fig_path}/mutation_probs_histogram_lin_scale.png")
        plt.clf()
        plt.close()


@click.command("examine_pcm_prob")
@click.option("-d", "--data-path", help="Path to data directory.")
# @click.option("-i", "--input-name",  help="Path to sample file to concatenate after clade, trial, sim_type.")
def examine_pcm_prob(data_path):
    """
    For creating plots that help us understand the relationship between a mutation's mutation rate
    and its probability of being involved in a PCM.
    """
    file_name = "counts_rate_list"  # Name of the file that we want to use for our scatter plot
    clade_names = [
        "AY.34.2",
        "AY.108",
        "AZ.3",
        "P.1.7",
        "AY.74",
        "A.2.5",
        "AY.87",
        "B.1.1.10"
    ]

    for clade_name in clade_names:
        clade_data = []
        max_mutation_rate = 0
        clade_x, clade_y = [], []
        clade_fig_path = f"{data_path}/{clade_name}/gamma_10_hmut_50/figures/diff_historydag_mut_rates"
        if not os.path.isdir(clade_fig_path):
                os.makedirs(clade_fig_path)

        for trial in range(1, 26):
            if trial == 1:
                print(clade_name)
            base_path = f"{data_path}/{clade_name}/gamma_10_hmut_50/{trial}"

            scatter_list_path = base_path + f"/results/diff_historydag_mut_rates/{file_name}.pkl"

            try:
                with open(scatter_list_path, "rb") as f:
                    scatter_list = pickle.load(f)
            except Exception as e:
                print(e)
                continue

            # Gather prob pcm info
            y = [num_pcm / count for num_pcm, count, mut_rate in scatter_list]
            x = [mut_rate for num_pcm, count, mut_rate in scatter_list]
            clade_x.extend(x)
            clade_y.extend(y)
            clade_data.extend(scatter_list)
            max_mutation_rate = max(max_mutation_rate, max(x))

            # NOTE: Can plot individual trees by uncommenting
            # trial_fig_path = base_path + f"/figures/diff_historydag_mut_rates"
            # if not os.path.isdir(trial_fig_path):
            #     os.makedirs(trial_fig_path)
            # plt.scatter(x, y)
            # plt.xlabel("Number of PCMs")
            # plt.ylabel("Mutation Rate")
            # plt.savefig(f"{trial_fig_path}/num_pcm_vs_rate.png")
            # plt.clf()
            # plt.close()

        # Plot the proportion of PCMs as a function of mutation rate
        plt.scatter(clade_x, clade_y, alpha=0.5)
        plt.ylabel("Number of PCMs")
        plt.xlabel("Mutation Rate")
        plt.savefig(f"{clade_fig_path}/{file_name}_clade_scatter.png")
        plt.clf()
        plt.close()


        num_bins = 4
        # NOTE: Could also try using geometrically-spaced bins
        # mr_bins = [i * max_mutation_rate / num_bins for i in range(num_bins)]
        # mr_bins.append(max_mutation_rate)
        mr_bins = [max_mutation_rate / 2 ** i for i in range(num_bins)]
        mr_bins.append(0)
        mr_bins.reverse()
        binned_data = {i: [] for i in range(num_bins)}
        for num_pcm, count, mut_rate in clade_data:
            filter_arr = [i for i, mr_bin in enumerate(mr_bins) if mr_bin < mut_rate]
            bin_number = filter_arr[-1]
            binned_data[bin_number].append((num_pcm, count))
        
        # Plot the mutation rate bins on different panels
        fig, axes = plt.subplots(num_bins, sharex=True)
        fig.set_size_inches(6, 10)
        fig.suptitle('PCM Counts Stratified by Mutation Rate')
        for i, ax in enumerate(axes):
            idx = num_bins-1-i
            data = binned_data[idx]
            x = [num_pcm for num_pcm, count in data]
            y = [count for num_pcm, count in data]
            ax.scatter(x, y, alpha=0.5)
            ax.set_title(f"MR in [{mr_bins[idx]:.2f}, {mr_bins[idx+1]:.2f}]")
            ax.set(ylabel='Count observed')
        axes[-1].set(xlabel='# times in a PCM')
        plt.savefig(f"{clade_fig_path}/{file_name}_mut_rate_stratified_scatter.png")
        plt.clf()
        plt.close()

        # Same panel with different colors
        plt.title('PCM Counts Stratified by Mutation Rate')
        for i in range(num_bins):
            data = binned_data[i]
            x = [num_pcm for num_pcm, count in data]
            y = [count for num_pcm, count in data]
            plt.scatter(x, y, alpha=0.4, c=colors[i], label=f"[{mr_bins[i]:.2f}, {mr_bins[i+1]:.2f}]")
        plt.xlabel('# times in a PCM')
        plt.ylabel('Count observed')
        plt.legend()
        plt.savefig(f"{clade_fig_path}/{file_name}_mut_rate_color_stratified_scatter.png")
        plt.clf()
        plt.close()




#####
#####
cli.add_command(examine_pcm_prob)
cli.add_command(compare_mutation_probs)
cli.add_command(plot_sub_stats)

cli.add_command(plot_pars_distribution)
cli.add_command(plot_sim_stats)   # NOTE: This should probably go in the hmut stuff...
cli.add_command(sim_model_coverage)
cli.add_command(plot_sim_model_stats)

cli.add_command(hmut_coverage)
cli.add_command(plot_hmut)

cli.add_command(coverage_trial_plot)
cli.add_command(compare_trial_results)
cli.add_command(compare_results)
cli.add_command(clade_results)
cli.add_command(compare_clade_results)
cli.add_command(aggregate_results)


###################################################################################################
### NOTE: Methods below here have not been thoroughly checked
###################################################################################################

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
    results_full = get_results_general(clade_dir, num_sim, method, results_name)
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

    

if __name__ == '__main__':
    cli()
