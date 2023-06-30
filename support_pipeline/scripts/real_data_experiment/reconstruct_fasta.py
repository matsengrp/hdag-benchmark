import bte
import sys


def reconstruct_fasta(clade_path):
    """
    Reconstructs 'original' fasta file from tree and leaf sequences. Expects USHER protobuf to
    be stored at `clade_tree.pb.gz` and the reference sequence to be at `refseq.fasta` both in the
    `clade_path` directory. Also, assumes that tree.mapping is in clade_path directory.

    NEED to run `conda activate bte` first
    """

    mapping = get_mapping(f"{clade_path}/tree.mapping")

    import bte
    tree = bte.MATree(f'{clade_path}/clade_tree.pb.gz')

    # Build new reference sequence
    with open(f'{clade_path}/../refseq.fasta', 'r') as f:
        f.readline()
        refseq_orig = f.readline()
    mutations = {}
    for node in tree.breadth_first_expansion():
        for mut in node.mutations:
            base = mut[0]
            idx = int(mut[1:-1])-1
            if idx not in mutations:
                mutations[idx] = base
    refseq = edit_string(refseq_orig, mutations)

    # Create edit dictionary for each node and apply edits
    leaf2seq = {}
    leaves = tree.get_leaves()
    for leaf in leaves:
        mutations = []
        mutations.extend(leaf.mutations)
        node = leaf
        while node.level > 1:
            node = node.parent
            mutations.extend(node.mutations)
            if node.level == 1:
                break
        
        mutations.reverse()
        node_seq = refseq
        muts = {}
        for mut in mutations:
            idx = int(mut[1:-1])-1
            new_base = mut[-1]
            muts[idx] = new_base
        node_seq = edit_string(refseq, muts)
        
        # Include each sequence once
        if node_seq not in leaf2seq.values():
            leaf2seq[leaf.id] = node_seq
    
    out_path = f"{clade_path}/reconstructed_seqs.fasta"
    with open(out_path, 'w') as f:
        f.write(f">ancestral\n{refseq_orig}")
        for leaf_id, seq in leaf2seq.items():
            f.write(f">{mapping[leaf_id]}\n{seq.strip()}\n")

def edit_string(original, mutations):
    """
    Given the original string and a list of mutations, returns the new edited string.
    """
    new_str = ""
    prev_idx = 0
    idx_list = list(mutations.keys())
    idx_list.sort()
    for idx in idx_list:
        new_str += original[prev_idx:idx] + mutations[idx]
        prev_idx = idx+1
    if idx+1 < len(original):
        new_str += original[prev_idx:len(original)+1]
    
    # Check to make sure your new sequence is correct
    assert len(new_str) == len(original)
    for i in range(len(original)):
        if i not in mutations.keys():
            assert original[i] == new_str[i]
        else:
            new_str[i] == mutations[i]

    return new_str

def get_mapping(fp, delim=" "):
    """
    Assumes file has lines of space separated key-value pairs.
    """
    with open(fp, "r") as f:
        mapping = {}
        for line in f.readlines():
            pair = line.strip().split(delim)
            mapping[pair[0]] = pair[1]
        return mapping


def main():
    reconstruct_fasta(sys.argv[1])

if __name__ == '__main__':
    main()

            
                