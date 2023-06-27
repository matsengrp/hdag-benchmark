
import ete3 as ete
import historydag as hdag
import time
import sys
import random

seed = sys.argv[1]
random.seed(seed)

base_dir = f"/fh/fast/matsen_e/whowards/hdag-benchmark/data_new_sim/AY.108/real"
out_path = base_dir + f"/{seed}/collapsed_hdag_tree.nwk"
dag_path = base_dir + "/results/historydag/final_opt_dag.pb"
# dag_path = "/fh/fast/matsen_e/whowards/hdag-benchmark/data_new_sim/AY.108/real/results/historydag/opt_dag_1.pb"

fasta_path = base_dir + "/collapsed_hdag_tree.nwk.fasta"


def main():
    start = time.time()
    print(f"Loading DAG from {dag_path}...")
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file(dag_path)
    print(f"\tLoaded in {time.time() - start} seconds")

    start = time.time()
    print("Trimming DAG...")
    dag.convert_to_collapsed()
    dag.make_complete()
    dag.convert_to_collapsed()
    dag.trim_optimal_weight()
    print(f"\tTrimmed in {time.time() - start} seconds")

    dag.summary()

    start = time.time()
    print("Sampling tree...")
    dag.uniform_distribution_annotate()
    sample = dag.sample()

    id2seq = hdag.parsimony.load_fasta(fasta_path)
    seq2id = {v:k for k,v in id2seq.items()}
    
    name_func = get_name_mapping(sample, seq2id)
    ete_tree = sample.to_ete(name_func = lambda n: name_func[n] if n.is_leaf() else 'unnamed')
    ete_tree.write(outfile=out_path, format=9)
    print(f"\tSampled in {time.time() - start} seconds")


def get_name_mapping(tree, seq2id):
    """
    Given sequence to name dictionary, and history from MADAG, returns a map of histordag leaf
    nodes to their name in the given fasta file.
    """

    node2name = {}
    for node in tree.postorder():
        if node.is_leaf():
            cg = node.label.compact_genome
            seq = cg.to_sequence()
            name = seq2id[seq]
        else:
            name = 'unnamed'
        node2name[node] = name

    return node2name


if __name__ == "__main__":
    main()