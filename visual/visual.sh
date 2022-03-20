#!/usr/bin/env bash
# Takes a vcf or fasta file as input 1 and output dir as optional input 2.
# Output is a figure showing the gene map of the input.
# Plots are creating with R scripts from https://github.com/jonas-fuchs/SARS-CoV-2-analyses .

input=$1

# check if outdir was provided and set if possible
if [[ $2 ]]; then
    output=$2
    mkdir -p $output
else
    output="../ouput"
    mkdir $output
fi

# check if input is vcf or fasta
if [ "${input##*.}" = "vcf" ]; then
    # script expects vcf files to be in same dir
    cp $input ../bin
    Rscript ../bin/Genome_variant_vis.R
    rm ../bin/*vcf
    rm Genome_variant_vis.pdf
    mv *.pdf $output
else
    # check if covsonar is available
    if ! python3 sonar.py -h; then
        git clone https://gitlab.com/s.fuchs/covsonar
        conda env create -n sonar -f covsonar/sonar.env.yml
        conda activate sonar
    fi

fi


#conda deactivate