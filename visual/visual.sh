#!/usr/bin/env bash
# Takes a vcf or fasta file as input 1 and output dir as optional input 2.
# Output is a figure showing the gene map of the input.
# Plots are creating with R scripts from https://github.com/jonas-fuchs/SARS-CoV-2-analyses .

input=$1
filename=${input##*/}
basename=${filename%.*}

# check if outdir was provided and set if possible
if [[ $2 ]]; then
    output=$2
else
    output="../ouput"
fi

    mkdir -p $output/vcf
    mkdir -p $output/figures

# check if input is a vcf file
if [ "${input##*.}" = "vcf" ]; then
    # script expects vcf files to be in same dir
    cp $input ../bin
    echo "running Rscript.."
    Rscript ../bin/Genome_variant_vis.R
    rm ../bin/*vcf
    #rm Genome_variant_vis.pdf
    mv *.pdf $output/figures
# assume input is a fasta file
else
    CONDA_BASE=$(conda info --base)
    # check if covsonar is available and set it up if not
    if [ ! -d covsonar ]; then
        echo "Covsonar is not available but a fasta file was provided as input."
        git clone https://gitlab.com/s.fuchs/covsonar
    fi
    if conda env list | grep sonar > /dev/null 2>&1; then
        source $CONDA_BASE/etc/profile.d/conda.sh
        conda activate sonar
        echo "current environment: $CONDA_DEFAULT_ENV"
    else
        source $CONDA_BASE/etc/profile.d/conda.sh
        conda env create -n sonar -f covsonar/sonar.env.yml
        conda activate sonar
        echo "current environment: $CONDA_DEFAULT_ENV"
    fi

    echo "running covsonar.."
    # adding all sequences from input to database 'mydb'
    # using eight cpus (the database file will be created if it does not exist)
    python3 covsonar/sonar.py add -f $input --db mydb --cpus 8
    # get vcf file from covsonar db
    python3 covsonar/sonar.py var2vcf --db mydb -o $output/vcf/merge.vcf --cpus 8
    rm -f mydb
    gunzip $output/vcf/merge.vcf.gz
    mv $output/vcf/merge.vcf $output/vcf/$basename.vcf

    echo "running Rscript.."
    cp $output/vcf/$basename.vcf ../bin
    Rscript ../bin/Genome_variant_vis.R
    rm ../bin/*vcf
    #rm Genome_variant_vis.pdf
    mv *.pdf $output/figures
fi

# conda env is deactivated when this sub-shell dies
#conda deactivate