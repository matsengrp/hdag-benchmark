#!/bin/bash
# Shell script that uses dnapars to compute the best parsimony score for
#   a tree constructed on a fasta.
# usage: bash dnapars_parsimony_score.sh input.fasta root_sequence_name output_directory

# provide input fasta, name of root sequence, and output directory to script

set -eu
infasta=$1
rootname=$2
outdir=$3
prefix=$outdir/dnapars

rm -rf $prefix
mkdir -p $prefix

# TODO: Test that outdir works

deduplicate $infasta --root $rootname --idmapfile $prefix/idmap.txt > $prefix/deduplicated.phylip
realpath $prefix/deduplicated.phylip > $prefix/dnapars.cfg
echo "J
1
10
O
1
4
5
.
V
1
Y
" >> $prefix/dnapars.cfg
old_dir=$(pwd)
cd $prefix
dnapars < dnapars.cfg > dnapars.log
cd $old_dir
grep "requires a total of" $prefix/outfile > $outdir/dnapars_output.txt