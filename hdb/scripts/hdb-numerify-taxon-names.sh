#!/bin/sh
#
# Convert taxon names to s<number>, outputting a "mapping file" and a renamed newick.

set -eu

test $# -ne 2 && {
    echo "Usage: hdb-numerify-taxon-names.sh tree.nwk"
    exit 1
}

orig_path=$1
out_path_base=$2
mapping_path=${out_path_base}.mapping
numerified=${out_path_base}.n.nwk

echo "ancestral ancestral" > $mapping_path
nw_labels $orig_path | grep -v ancestral | awk '{print $0 " s" NR}' >> $mapping_path
nw_rename $orig_path $mapping_path > $numerified
