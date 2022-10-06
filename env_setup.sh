set -eu

cd data
wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2022/10/01/public-2022-10-01.all.masked.pb.gz
python choose_clades.py > clades.txt
for clade in $(cat clades.txt); do
    # send a cluster job off that does everything to that clade?

