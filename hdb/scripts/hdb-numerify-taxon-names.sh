#!/bin/sh
#
# Convert taxon names to s<number>, outputting a "mapping file" and a renamed newick.

set -eu

test $# -ne 1 && {
    echo "Usage: hdb-numerify-taxon-names.sh tree.nwk"
    exit 1
}

orig_path=$1
base=$(basename $orig_path .nwk)
mapping_path=${base}.mapping
numerified=${base}.n.nwk

echo "ancestral ancestral" > $mapping_path
nw_labels $orig_path | grep -v ancestral | awk '{print $0 " s" NR}' >> $mapping_path
nw_rename $orig_path $mapping_path > $numerified
