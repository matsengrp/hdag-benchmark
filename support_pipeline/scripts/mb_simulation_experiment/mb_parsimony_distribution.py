import historydag as hdag
from historydag.parsimony import load_fasta, build_tree, disambiguate, sankoff_upward
import re
import sys
import pickle

experiment_type = "sim_models"
sim_type = sys.argv[1]
clade = sys.argv[2]
trial = sys.argv[3]


# Number of trees in the MrBayes run
num_trees = 5e4
burnin = 0.9
iter_per_sample = 1000

dir_path = f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/{experiment_type}/{clade}/{sim_type}/{trial}"

fasta_path = dir_path + "/simulation/collapsed_simulated_tree.nwk.fasta"
mb_path = dir_path + "/results/mrbayes/mrbayes-output.t"
fasta = load_fasta(fasta_path)

def main():
    distr = get_parsimony_distribution(mb_path)
    print(distr)
    with open(dir_path + "/results/mrbayes/mrbayes-pars-distr.pkl", "wb") as f:
        pickle.dump(distr, f)


def get_trees(input_path):
    i = -1
    leaf_map = {}
    with open(input_path, 'r') as fh:
        for line in fh:
            if '    ' in line[0:5]:
                mb_leaf, fasta_leaf = tuple(line.strip().split(' '))
                leaf_map[mb_leaf] = fasta_leaf[:-1]
                fl2ml = {v: k for k, v in leaf_map.items()}
            elif 'tree' in line.strip()[0:5]:
                i += 1
                if i < burnin * num_trees:
                    continue
                treeprob = 1
                nwk = line.strip().split(' ')[-1]

                fasta_mb = {fl2ml[fasta_leaf]: v for fasta_leaf, v in fasta.items()}
                tree = build_tree(nwk, fasta_mb)

                pars_score = sankoff_upward(
                    tree,
                    len(tree.sequence),
                    transition_model=hdag.parsimony_utils.default_nt_transitions,
                    use_internal_node_sequences=False,
                )
                tree.pars_score = pars_score
                tree.prob = treeprob

                if i % 10 == 0:
                    print(i, pars_score)

                yield tree
                
def get_parsimony_distribution(input_path):
    distribution = dict()
    for tree in get_trees(input_path):
        pars = tree.pars_score
        if pars in distribution:
            distribution[pars] += tree.prob
        else:
            distribution[pars] = tree.prob
    return distribution


if __name__ == '__main__':
    main()