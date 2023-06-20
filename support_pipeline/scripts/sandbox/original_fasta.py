import bte


def reconstruct_fasta():
    """
    Reconstructs 'original' fasta file from tree and leaf sequences

    NEED to run `conda activate bte` first
    """
    tree = bte.MATree('/fh/fast/matsen_e/whowards/hdag-benchmark/data_new_sim/AY.108/real/clade.pb.gz')

    # Build new reference sequence
    with open('/fh/fast/matsen_e/whowards/hdag-benchmark/data_new_sim/refseq.fasta', 'r') as f:
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
    
    out_path = "/fh/fast/matsen_e/whowards/hdag-benchmark/data_new_sim/AY.108/real/real_seqs.fasta"
    leaf2seq["ancestral"] = refseq_orig
    with open(out_path, 'w') as f:
        for leaf_id, seq in leaf2seq.items():
            f.write(f">{leaf_id}\n{seq.strip()}\n")

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

def main():
    reconstruct_fasta()

if __name__ == '__main__':
    main()

            
                