import historydag as hdag
from pathlib import Path
from historydag.parsimony import load_fasta, build_tree, disambiguate, sankoff_upward
import time
import sys
import pickle
import os
import ete3

experiment_type = "sim_models"
sim_type = "gamma_10_hmut_50"
clade = sys.argv[1]
trial = sys.argv[2]
rf_method = sys.argv[3]

input_dir_path = f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/{experiment_type}/{clade}/{sim_type}/{trial}"
out_dir_path = f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/rf_scaling/{clade}/{trial}/results"

def main():
    stats = get_stats(input_dir_path)
    print(stats)
    os.makedirs(out_dir_path, exist_ok=True)
    with open(out_dir_path + f"/{rf_method}.pkl", "wb") as f:
        pickle.dump(stats, f)

                
def get_stats(input_dir_path):
    """
    Given a path to a hdag protobuf, returns a dictionary containing statistics about the number
    of histories in the hdag and the time it takes to run naive inference
    """
    stats = {
        "number_histories": None,
        "clock_time": None,
        "number_leaves": None
    }

    hdag, sim_history = load_dags(input_dir_path)

    stats["number_histories"] = hdag.count_histories()
    stats["number_leaves"] = len(list(hdag.get_leaves()))

    start = time.time()
    if rf_method == "naive":
        tree = naive_min_tree(hdag, sim_history)
    elif rf_method == "min_trim":
        tree = hdag_min_tree(hdag, sim_history)
    else:
        raise Warning("Invalid minimizer method!")
    end = time.time()
    stats["clock_time"] = end - start

    return stats

def naive_min_tree(hdag, sim_history):
    min_dist_history = None
    min_dist = 1e20
    for history in hdag.get_histories():
        dist = history.optimal_rf_distance(sim_history)
        if dist < min_dist:
            min_dist = dist
            min_dist_history = history
    return min_dist_history

def hdag_min_tree(hdag, sim_history):
   hdag.trim_optimal_rf_distance(sim_history)
   return hdag.sample()

def load_dags(input_dir_path):
    """
    Returns the collapsed simulated history and the historydag we build with larch for this clade-trial
    """
    inpath = Path(input_dir_path)

    print("Input path:", str(inpath))

    dagpath = inpath / 'results/historydag/final_opt_dag_trimmed.pb'
    sim_inpath = inpath / 'simulation'
    fasta = hdag.parsimony.load_fasta(sim_inpath / 'ctree_with_refseq.fasta')
    etetree = ete3.Tree(str(sim_inpath / 'collapsed_simulated_tree.nwk'), format=0)

    # Reconstruct sequences on each node
    ancestral_seq = fasta["ancestral"]
    for node in etetree.traverse("preorder"):
        if node.is_root():
            assert len(node.mutations) == 0
            parent_seq = ancestral_seq
        else:
            parent_seq = node.up.sequence
        
        if len(node.mutations) > 0:
            mut_dict = {}
            mut_list = node.mutations.split("|")
            for mut_str in mut_list:
                parent_seq = replace_site(int(mut_str[1:-1]), mut_str[-1], parent_seq, oldbase=mut_str[0])
    
        node.add_feature("sequence", parent_seq)

    seqs = set(fasta.values())
    # NOTE: Debug check.
    for leaf in etetree.get_leaves():
        if leaf.sequence not in seqs:
            print("Mutations on naughty leaf:")
            mut_dict = {}
            curr = leaf
            while curr.up is not None:
                print("\t", curr.mutations.strip().split("|"))
                if len(curr.mutations) > 0:
                    for mut in curr.mutations.strip().split("|"):
                        idx = int(mut[1:-1])
                        from_nuc = mut[0]
                        to_nuc = mut[-1]
                        if idx not in mut_dict:
                            mut_dict[idx] = []
                        mut_dict[idx].append((from_nuc, to_nuc))
                curr = curr.up
            
            for idx, pair_list in mut_dict.items():
                print(f"{idx}\t{pair_list}")

    sim_history = hdag.history_dag_from_trees(
        [etetree],
        [],
        label_functions = {"sequence": lambda n: n.sequence},
        attr_func = lambda n: {"name": n.name if n.is_leaf() else "internal"},
    )
    sim_history.convert_to_collapsed()

    opt_dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(dagpath)
    opt_dag.convert_to_collapsed()

    for leaf in opt_dag.get_leaves():
        seq = leaf.label.compact_genome.to_sequence()
        if seq not in fasta.values():
            raise Warning("DAG leaves are not the same as fasta leaves. Aborting.")
        
    return opt_dag, sim_history

def replace_site(site, newbase, sequence, oldbase=None):
    idx = site - 1
    if oldbase is not None and sequence[idx] != oldbase:
        print(f"WARNING: previous base was not as expected (expected={oldbase} vs actual={sequence[idx]} at {site})")
    return sequence[:idx] + newbase + sequence[idx + 1:]

if __name__ == '__main__':
    main()