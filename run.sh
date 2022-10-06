tree=example_tree.nwk
cd data
mkdir _ignore
hdb-numerify-taxon-names.sh $tree
phastSim --outpath _ignore/ --seed 9 --createFasta --createInfo \
    --createNewick --createPhylip --treeFile $(basename -s .nwk $tree).n.nwk \
         --scale 0.00004 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 \
         --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon \
         --reference ref.fasta
hdb alnsummarize -o _ignore/summary.csv _ignore/sars-cov-2_simulation_output.fasta ../data/unique_seqs.fasta
column -t -s, _ignore/summary.csv| less
cd ..
