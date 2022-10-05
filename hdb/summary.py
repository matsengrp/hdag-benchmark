"""Alignment summary statistics."""

import pandas as pd
from amas import AMAS

def summary_df_of_fasta_paths(fasta_paths):
    "Get an alignment summary via AMAS."
    meta_aln = AMAS.MetaAlignment(in_files=fasta_paths, data_type="dna",in_format="fasta", cores=1)
    (headers, summaries) = meta_aln.get_summaries()
    df = pd.DataFrame(summaries, columns=headers)
    df.set_index("Alignment_name", inplace=True)
    return df
