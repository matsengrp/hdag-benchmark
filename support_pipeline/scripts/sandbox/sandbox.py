import historydag as hdag
from historydag.parsimony import load_fasta, build_tree, disambiguate, sankoff_upward
import re
import sys
import pickle

trial = sys.argv[1]

dir_path = f"/fh/fast/matsen_e/whowards/hdag-benchmark/data/AY.108/{trial}"
fasta_path = dir_path + "/simulation/collapsed_simulated_tree.nwk.variant_sites_with_refseq.fasta"
mb_path = dir_path + "/results/inference__mrbayes.log" # TODO: Change this to `inference__mrbayes.log`` for updated version

fasta = load_fasta(fasta_path)


def main():
    distr = get_distr(mb_path)
    min_pars = min(list(distr.keys()))
    max_pars = max(list(distr.keys()))
    cdf = [(min_pars, distr[min_pars])]
    for i in range(min_pars+1, max_pars+1):
        cdf.append((i, distr[i]+cdf[-1][1]))

    for p, score in cdf:
        print(p, "\t", score)


def get_distr(input_path):
    counts = {}
    with open(input_path, 'r') as fh:
        for line in fh:
            if '{' in line[0:1]:
                line = line[1:-2].split(', ')
                total = 0
                for txt in line:
                    k = int(float(txt.split(': ')[0].strip()))
                    v = int(txt.split(': ')[1].strip())
                    counts[k] = v
                    total += v
                distr = {pars: count / total for pars, count in counts.items()}
                return distr



if __name__ == '__main__':
    main()