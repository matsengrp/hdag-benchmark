import random
import pickle
from collections import Counter
from math import isclose
import numpy as np

import matplotlib.pyplot as plt
from scipy.stats import linregress, t

def sliding_window_plot(results, std_dev=False, sup_range=False, window_size=200, mean_x=False):
    """Given list of results tuples returns xy coords of sliding window plot."""

    results.sort(key= lambda result: result[1])

    x, y = [], []
    devs, min_sup, max_sup = [], [], []
    side_len = int(window_size/2)
    for i, (_, est_sup, in_tree) in enumerate(results):
        if i % (len(results)//10) == 0:
            print("\t", i)

        if mean_x:
            x_window = [float(el[1]) for el in results[max(0, i-side_len):min(len(results), i+side_len)]]
            x.append(sum(x_window) / len(x_window))
        else:
            x.append(est_sup)

        window = [int(el[2]) for el in results[max(0, i-side_len):min(len(results), i+side_len)]]
        y.append(sum(window) / len(window))
        
        # Standard deviations are over the window of in_true_tree variable (ie 1 or 0)
        if std_dev:
            avg = y[i]
            sq_diff = [(el - avg)**2 for el in window]
            devs.append((sum(sq_diff) / len(sq_diff)))

        # Quartiles are between estimated support values
        elif sup_range:
            support_window = [el[1] for el in results[max(0, i-side_len):min(len(results), i+side_len)]]
            # First and fourth quartiles
            if i - len(support_window)/2 <= 0:
                min_sup.append(support_window[0])
            else:
                min_sup.append(support_window[int(len(support_window)/4)])

            if i + len(support_window)/2 >= window_size:
                max_sup.append(support_window[-1])
            else:
                max_sup.append(support_window[-int(len(support_window)/4)])

    if std_dev:
        pos_devs = [y_val + dev for y_val, dev in zip(y, devs)]
        neg_devs = [y_val - dev for y_val, dev in zip(y, devs)]
        return x, y, pos_devs, neg_devs
    elif sup_range:
        return x, y, min_sup, max_sup
    else:
        return x, y


# def sliding_window_plot_new(results, std_dev=False, sup_range=False, window_size=200):
#     """Given list of results tuples returns xy coords of sliding window plot."""

#     results.sort(key= lambda result: result[1])

#     x, y = [], []
#     devs, min_sup, max_sup = [], [], []
#     side_len = int(window_size/2)
#     true_vals = Counter(res[2] for res in results[0 : min(len(results), side_len)])
#     sum_estimates = sum(res[1] for res in results[0: min(len(results), side_len)])
#     first_idx = 0
#     first_full_window_end_idx = min(len(results), window_size - 1)

#     for central_idx in range(len(results)):
#         proposed_first_idx = central_idx - side_len
#         proposed_last_idx = central_idx + side_len

#         # Adjust the window appropriately
#         if proposed_first_idx > first_idx:
#             assert proposed_first_idx - first_idx == 1
#             sum_estimates -= results[first_idx][1]
#             true_vals.subtract([results[first_idx][2]])
#             first_idx = proposed_first_idx
#         if proposed_last_idx < len(results):
#             sum_estimates += results[proposed_last_idx][1]
#             true_vals.update([results[proposed_last_idx][2]])
#         last_idx = min(proposed_last_idx, len(results))

#         this_window_size = sum(true_vals.values())
#         x.append(results[central_idx][1])
#         # #This would do averages instead of medians, but it doesn't seem to
#         # #work correctly..
#         # x.append(sum_estimates / this_window_size)
#         y.append(true_vals[1] / this_window_size)
#         # Standard deviations are over the window of in_true_tree variable (ie 1 or 0)
#         if std_dev:
#             avg = y[-1]
#             sq_diff = [el_count * (el - avg)**2 for el, el_count in true_vals.values()]
#             devs.append((sum(sq_diff) / this_window_size))

#         # Quartiles are between estimated support values
#         elif sup_range:
#             # First and fourth quartiles
#             min_sup.append(results[first_idx + int(this_window_size/4)])
#             max_sup.append(support_window[last_idx - int(this_window_size/4)])

#     if std_dev:
#         pos_devs = [y_val + dev for y_val, dev in zip(y, devs)]
#         neg_devs = [y_val - dev for y_val, dev in zip(y, devs)]
#         return x, y, pos_devs, neg_devs
#     elif sup_range:
#         return x, y, min_sup, max_sup
#     else:
#         return x, y

# def sliding_window_plot_new?(results, std_dev=False, sup_range=False, window_size=200):
#     """Given list of results tuples returns xy coords of sliding window plot."""
#     side_len = window_size // 2
#     window_size = side_len * 2 + 1
#     if len(results) < window_size:
#         raise ValueError("too few clades to compute sliding window plot")

#     results.sort(key= lambda result: result[1])

#     x, y = [], []
#     devs, min_sup, max_sup = [], [], []
#     true_vals = Counter(res[2] for res in results[0 : window_size])
#     sum_estimates = sum(res[1] for res in results[0: window_size])
#     center_idx = side_len
#     # These are the indices of the first and last records in the window
#     # (so the window is results[first_idx:last_idx + 1])
#     first_idx = 0
#     last_idx = window_size - 1

#     def log_stats(true_vals, sum_estimates, center_idx, first_idx, last_idx):
#         x.append(results[center_idx][1])
#         y.append(true_vals[1] / window_size)

#         observed_window = [it[2] for it in results[first_idx: last_idx + 1]]
#         estimates_window = [it[1] for it in results[first_idx: last_idx + 1]]
#         assert last_idx - first_idx == window_size - 1, f"{last_idx - first_idx} != {window_size - 1} (expected)"
#         assert isclose(y[-1], sum(observed_window) / window_size)


#     while last_idx < len(results) - 1:
#         log_stats(true_vals, sum_estimates, center_idx, first_idx, last_idx)
#         last_idx += 1
#         true_vals.update([results[last_idx][2]])
#         true_vals.subtract([results[first_idx][2]])
#         sum_estimates += results[last_idx][1] - results[first_idx][1]
#         first_idx += 1
#         center_idx += 1
#     log_stats(true_vals, sum_estimates, center_idx, first_idx, last_idx)


#     return x, y
#     # #### past here are old things that need to get integrated
#     #     if std_dev:
#     #         avg = y[-1]
#     #         sq_diff = [el_count * (el - avg)**2 for el, el_count in
#     #         true_vals.items()]
#     #         devs.append((sum(sq_diff) / this_window_size))

#     #     # Quartiles are between estimated support values
#     #     elif sup_range:
#     #         # First and fourth quartiles
#     #         min_sup.append(results[first_idx + int(this_window_size/4)])
#     #         max_sup.append(support_window[last_idx - int(this_window_size/4)])

#     # if std_dev:
#     #     pos_devs = [y_val + dev for y_val, dev in zip(y, devs)]
#     #     neg_devs = [y_val - dev for y_val, dev in zip(y, devs)]
#     #     return x, y, pos_devs, neg_devs
#     # elif sup_range:
#     #     return x, y, min_sup, max_sup
#     # else:
#     #     return x, y


def plot_scatter_with_line(df, x_col, y_col, save_path):
    # Extract the data from the dataframe
    x = df[x_col]
    y = df[y_col]

    idxs = np.argsort(x)
    x = x[idxs]
    y = y[idxs]

    # Perform Pearson correlation test
    correlation = np.corrcoef(x, y)[0, 1]
    slope, intercept, r, p_value, se = linregress(x, y)

    # Calculate the line of best fit
    line_of_best_fit = slope * x + intercept

    # Calculate confidence intervals
    y_err = y - line_of_best_fit
    mean_x = np.mean(x)
    n = len(x)
    dof = n - 2  # degrees of freedom
    alpha = 0.05  # Significance level (1 - confidence level)
    t_critical = t.ppf(1 - alpha / 2, dof)  # Calculate the critical value
    confidence_interval = t_critical * np.sqrt((np.sum(y_err ** 2) / dof) * (1.0 / n + (x - mean_x) ** 2 / np.sum((x - mean_x) ** 2)))

    # Plot the scatter plot and line of best fit
    plt.clf()
    plt.scatter(x, y, label='Data Points', alpha=0.1)
    plt.plot(x, line_of_best_fit, color='red', label='Line of Best Fit')

    # Plot confidence intervals
    plt.fill_between(x, line_of_best_fit - confidence_interval, line_of_best_fit + confidence_interval, color='orange', alpha=0.4, label='Confidence Intervals')

    # Set plot labels and title
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title('Scatter Plot with Line of Best Fit')

    # Add correlation coefficient, p-value, and critical value to the plot
    # text = f'Correlation: {correlation:.2f}\np-value: {p_value:.4f}\nCritical Value: {t_critical:.4f}'
    # plt.text(0.95, 0.05, text, transform=plt.gca().transAxes, ha='right', va='bottom', bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round'))

    # Add legend
    # plt.legend()

    # Save the plot
    plt.savefig(save_path)

    # Show the plot (optional)
    plt.show()

    
def get_results_full(clade_dir, num_sim, method, results_name, skip_list=[], support_removal_threshold=0.0):
    """
    Helper method for `clade_results` that aggregates all the results.pkl files for each trial
    into a single sorted list.

    Removes any node with estimated support that is less than or equal to `support_removal_threshold`.
    """
    raise Warning("Deprecating in favor of get_results_general")

    result_dict = {}
    for trial in range(1, num_sim+1):
        # Assumes that `path/to/clade/trial/results/results.pkl stores`` list of nodes
        #   their supports and whether they're in the true tree or not
        result_path = clade_dir + f"/{trial}/results/{method}/{results_name}"

        # Skip any trials you don't want (e.g., inference on them hasn't finished for this run)
        if trial in skip_list:
            continue
        
        try:
            with open(result_path, "rb") as f:
                results = pickle.load(f)

                # NOTE: Removing leaves and UA node here
                with_leaves = len(results)
                leaf_in_tree = [int(result[2]) for result in results if len(result[0]) <= 1]
                leaf_est_sup = [result[1] for result in results if len(result[0]) <= 1]
                # print(leaf_in_tree[0:10])
                # print(leaf_est_sup[0:10])
                results = [result for result in results if len(result[0]) > 1]
                without_leaves = len(results)
                if with_leaves != without_leaves:
                    print(f"==> Removed {with_leaves - without_leaves} leaves \
                        avg in_tree = {sum(leaf_in_tree) / len(leaf_in_tree)} \
                        avg est_sup = {sum(leaf_est_sup) / len(leaf_est_sup)}")
                
                result_dict[trial] = results
        except:
            print(f"\tSkipping {clade_dir} {trial}")
            continue
    
    if len(result_dict) == 0:
        print("\n==>No results to print :(\n")
        return

    results_full = []
    for trial, results in result_dict.items():
        for result in results:
            node = result[0]
            if len(node) > 1:       # Removing leaves
                if result[1] <= support_removal_threshold:
                    continue
                results_full.append((result[0], result[1], result[2]))
    
    print(f"\tsorting {len(results_full)} results...")
    random.shuffle(results_full)
    results_full.sort(key=lambda el: el[1])
    return results_full


def get_results_general(base_dir, num_sim, method, results_name, skip_list=[], support_removal_threshold=0.0):
    """
    Helper method that aggregates all the results.pkl files for each trial into a single sorted list.

    Removes any node with estimated support that is less than or equal to `support_removal_threshold`.
    """

    result_dict = {}
    for trial in range(1, num_sim+1):
        # Assumes that `path/to/clade/trial/results/results.pkl stores`` list of nodes
        #   their supports and whether they're in the true tree or not
        result_path = base_dir + f"/{trial}/results/{method}/{results_name}"

        # Skip any trials you don't want (e.g., inference on them hasn't finished for this run)
        if trial in skip_list:
            continue
        
        try:
            with open(result_path, "rb") as f:
                results = pickle.load(f)

                # NOTE: Removing leaves and UA node here
                # with_leaves = len(results)
                # leaf_in_tree = [int(result[2]) for result in results if len(result[0]) <= 1]
                # leaf_est_sup = [result[1] for result in results if len(result[0]) <= 1]
                # print(leaf_in_tree[0:10])
                # print(leaf_est_sup[0:10])
                results = [res for res in results if len(res[0]) > 1 and res[1] > support_removal_threshold]
                # without_leaves = len(results)
                # if with_leaves != without_leaves:
                #     print(f"==> Removed {with_leaves - without_leaves} leaves \
                #         avg in_tree = {sum(leaf_in_tree) / len(leaf_in_tree)} \
                #         avg est_sup = {sum(leaf_est_sup) / len(leaf_est_sup)}... These should both be 1.0")
                
                result_dict[trial] = results
        except:
            print(f"\tSkipping {base_dir} {trial}")
            continue
    
    if len(result_dict) == 0:
        print("\n==>No results to print :(\n")
        return

    results_full = []
    for trial, results in result_dict.items():
        results_full.extend(results)#[res for res in results if len(res[0]) > 1 and res[1] > support_removal_threshold])
        # node = result[0]
        # if len(node) > 1:       # Removing leaves
        #     if result[1] <= support_removal_threshold:
        #         continue
        #     results_full.append((result[0], result[1], result[2]))
    
    print(f"\tsorting {len(results_full)} results...")
    random.shuffle(results_full)
    results_full.sort(key=lambda el: el[1])
    return results_full


def bin_hist_plot(results, bin_size=0.05):
    """Given list of results tuples returns xy coords of sliding window plot."""

    results.sort(key= lambda result: result[1])
    x, y, devs = [], [], []
    
    ind_max = 0
    for bin in range(0, 100, int(bin_size * 100)): # NOTE: this is in percent
        est_max = bin / 100 + bin_size

        ind_min = ind_max   # Use previous maximum index
        for i, el in enumerate(results[ind_min:]):
            if el[1] >= est_max and est_max != 1:
                break
        ind_max = ind_min + i
        if i == 0:
            continue

        in_tree_window = [int(el[2]) for el in results[ind_min:ind_max]]
        avg = sum(in_tree_window) / len(in_tree_window)
        y.append(avg)
        est_sup_window = [el[1] for el in results[ind_min:ind_max]]
        x.append(sum(est_sup_window) / len(est_sup_window))

        # Standard deviations are over the window of in_true_tree variable (ie 1 or 0)
        sq_diff = [(el - avg)**2 for el in in_tree_window]
        devs.append((sum(sq_diff) / len(sq_diff)))

    pos_devs = [y_val + dev for y_val, dev in zip(y, devs)]
    neg_devs = [y_val - dev for y_val, dev in zip(y, devs)]

    return x, y, pos_devs, neg_devs
