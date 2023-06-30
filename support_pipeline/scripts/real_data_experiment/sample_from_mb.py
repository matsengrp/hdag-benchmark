# Given a path to a MrBayes generated tree path, and a prefix to an out p

import sys
sys.path.append("/fh/fast/matsen_e/whowards/hdag-benchmark/support_pipeline")

from utils import get_mrbayes_trees
import random
random.seed(123)

def main(tree_file_path, output_path):
    if tree_file_path[-1] == "/":
        raise(Warning("Will append additional '/'"))
    if output_path[-1] == "/":
        raise(Warning("Will append additional '/'"))
    sample_topology(tree_file_path, output_path)

def sample_topology(tree_file_path, output_path, burnin=0.7, num_iter=5e7):
    # sample_freq = 1000
    # num_trees = num_iter / sample_freq
    
    for i, tree in enumerate(get_mrbayes_trees(tree_file_path)):
        continue
    num_trees = i
    idx = random.randint(0, int((num_trees - burnin * num_trees)))

    tree_count = 0
    tree_sampled = False
    for i, tree in enumerate(get_mrbayes_trees(tree_file_path)):
        if i < burnin * num_trees:
            continue
        if tree_count >= idx:
            tree.write(format=9, outfile=output_path)
            print(f"Writing tree to {output_path}\n", tree.write(format=9))
            tree_sampled = True
            break
        tree_count += 1

    if not tree_sampled:
        raise(Warning(f"No tree sampled. Sample index was {idx}"))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])