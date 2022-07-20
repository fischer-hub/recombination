#!/usr/bin/env bash

# This script runs the steps for pre-screening new sequences for possible 
# recombination canidates. Steps include: 
# 1) alignin to the reference sequence used by sc2rf using nextalign
# 2) running sc2rf on the sequence alignment
# 3) filtering the sc2rf output for possible candidates
# 4) subsampling these candidate sequences from the original fasta file
# 5) checking whether sequences among these candidates are already identified
# as recombs by pangolin
#
# To run the script just start it like a normal bash script, all dependencies
# should be installed automatically. If this doesn't work please create a
# conda environment from the recomb_check.yaml and run the script inside of
# the environment. You can find the output files in the directory 
# "output_$current_date" or define your own output directory via the $out
# variable below. The same applies for the input file ($input). Input file can
# now also be provided via the first script argument.

input=""

if [ -n "$1" ]; then
    in=$1
    echo "input (arg 1): $in"
else
    in=$input
    echo "input: $in"
fi
out="output_$(date +%F)"

# check if env is set up and set up if its not
if conda env list | grep recomb_check > /dev/null 2>&1; then
    echo "environment is available, activating .."
    source $(conda info --base)/etc/profile.d/conda.sh
    [ $"CONDA_DEFAULT_ENV" != "recomb_check" ] && conda activate recomb_check
    echo "current environment: $CONDA_DEFAULT_ENV"
else
    echo "installing environment .."
    source $(conda info --base)/etc/profile.d/conda.sh
    conda env create -n recomb_check -f recomb_check.yaml
    conda activate recomb_check
    echo "current environment: $CONDA_DEFAULT_ENV"
fi

# copy new desh fasta
cp $in .

# align sequences with nextalign 
if ! [ -x "$(command -v nextalign)" ]; then
    echo 'nextalign not available, trying to install now ..' >&2
    conda install -c bioconda nextalign
fi

if ! [ -d "sc2rf" ]; then
    echo 'sc2rf not available, trying to install now ..' >&2
    git clone https://github.com/lenaschimmel/sc2rf.git
fi

echo 'running nextalign ..'
nextalign run --input-ref sc2rf/reference.fasta --output-fasta day.aligned.fasta day.fasta

# run sequences through sc2rf
echo 'running sc2rf ..'
cd sc2rf
python3 sc2rf.py --csvfile ../day.csv --parents 1-35 --breakpoints 1-2 \
                       --max-intermission-count 3 --max-intermission-length 1 \
                       --unique 1 --max-ambiguous 10000 --max-name-length 55 \
                       ../day.aligned.fasta > /dev/null 2>&1
cd ..

# filter out sequences that have any other donor lineage except Omicron and Delta and extract their IDs
awk -F ',' '$2~/Delta|Omicron/ && $3~/Delta|Omicron/{if($4 <= 3 && $5 >=1 && $5 <= 2) print $1}' day.csv > potenitalRecombinants.ids

# subsample all possible candidates from original fasta file
echo "running seqkit (potenitalRecombinants) .."
seqkit grep -f potenitalRecombinants.ids day.fasta > potenitalRecombinants.fasta 2> /dev/null

# find recomb candidates that already have been identified by pangolin
echo "running pangolin .."
pangolin potenitalRecombinants.fasta --skip-scorpio -t 4 > /dev/null 2>&1;
awk -F ',' '$2~/X/{print $1}' lineage_report.csv > pangolinDetectedRecombinants.ids


# filter day.csv for known recomb candidate seqs
echo "filtering day.csv for known candidates .."
head -n 1 day.csv > pangolinDetectedRecombinants.csv
while IFS="" read -r p || [ -n "$p" ]
do
   grep "$p" day.csv >> pangolinDetectedRecombinants.csv
done < pangolinDetectedRecombinants.ids

# filter day.csv for recomb candidate seqs
echo "filtering day.csv for new candidates .."
head -n 1 day.csv > potenitalRecombinants.csv
while IFS="" read -r p || [ -n "$p" ]
do
   grep "$p" day.csv >> potenitalRecombinants.csv
done < potenitalRecombinants.ids

# filter ids for unidentified ids
echo "filtering all candidates for unidentified candidates only .."
grep -v -f pangolinDetectedRecombinants.ids potenitalRecombinants.ids >> undetectedRecombinants.ids
grep -v -f pangolinDetectedRecombinants.ids day.csv >> undetectedRecombinants.csv

#subsample the initial fasta file for those sequences only
echo "running seqkit (pangolinDetectedRecombinants) .."
seqkit grep -f pangolinDetectedRecombinants.ids day.fasta > pangolinDetectedRecombinants.fasta 2> /dev/null

echo "running seqkit (undetectedRecombinants) .."
seqkit grep -v -f pangolinDetectedRecombinants.ids potenitalRecombinants.fasta > undetectedRecombinants.fasta 2> /dev/null

# get lines from pangolin lineage report for pangolinDetectedRecombinants and undetectedRecombinants
cat lineage_report.csv | head -n 1 > lineage_report_detected_tmp
grep -f pangolinDetectedRecombinants.ids lineage_report.csv >> lineage_report_detected_tmp

cat lineage_report.csv | head -n 1 > lineage_report_undetected_tmp
grep -f undetectedRecombinants.ids lineage_report.csv >> lineage_report_undetected_tmp


echo -e "\nfound $(cat potenitalRecombinants.ids | wc -l) total recombinant candidates"
echo "of which $(cat pangolinDetectedRecombinants.ids | wc -l) have already been indentified by pangolin"
echo "and $(cat undetectedRecombinants.ids | wc -l) are unidentified"

# clean up
echo -e "\ncleaning up .."
mkdir -p $out
rm day.fasta day.aligned.fasta lineage_report.csv
mv lineage_report_detected_tmp lineage_report_pangolinDetectedRecombinants.csv
mv lineage_report_undetected_tmp lineage_report_undetectedRecombinants.csv
mv day.csv $out
mv *.ids $out
mv *.fasta $out
mv *.csv $out

echo "done!"
