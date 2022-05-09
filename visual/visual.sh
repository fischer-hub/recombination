#!/usr/bin/env bash

### Description:
# Takes a vcf or fasta file as input 1 and output dir as optional input 2.
# Output is a figure showing the gene map of the input.

### Credits:
# Plots are created with R scripts from https://github.com/jonas-fuchs/SARS-CoV-2-analyses .
# VCF files from FASTA are created with covSonar from https://gitlab.com/s.fuchs/covsonar .



input=${1%/} # Remove possible trailing /

for filepath in $input/*; do

    filename=${filepath##*/}
    basename=${filename%.*}
    echo "currently processing: ${filepath}"
    # check if outdir was provided and set if possible
    if [[ $2 ]]; then
        output=$2
    else
        output="../ouput"
    fi

        mkdir -p $output/vcf
        mkdir -p $output/figures

    # check if input is a vcf file
    if [ "${filename##*.}" = "vcf" ]; then
        # script expects vcf files to be in same dir
        cp $filepath ../bin
        echo "running Genome_variant_vis.R .."
        Rscript ../bin/Genome_variant_vis.R
        rm ../bin/*vcf
        mv *.pdf $output/figures
    # assume input is a fasta file
    elif [ "${filename##*.}" = "gff" ]; then
        echo "this is a .gff file, skipping .."
        continue
    elif [[ ! $filename =~ "." ]]; then
        echo "unexpected file, skipping .."
        continue
    else
        CONDA_BASE=$(conda info --base)
        # check if covsonar is available and set it up if not
        if [ ! -d covsonar ]; then
            echo "Covsonar is not available but a fasta file was provided as input."
            git clone https://gitlab.com/s.fuchs/covsonar
        fi
        # check if sonar conda env is avalable, set up if not
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
        python3 covsonar/sonar.py add -f $filepath --db mydb --cpus 8
        echo "python3 covsonar/sonar.py add -f $filepath --db mydb --cpus 8"
        # get vcf file from covsonar db
        python3 covsonar/sonar.py var2vcf --db mydb -o $output/vcf/merge.vcf --cpus 8
        rm -f mydb
        gunzip $output/vcf/merge.vcf.gz
        mv $output/vcf/merge.vcf $output/vcf/$basename.vcf

        echo "running Genome_variant_vis.R .."
        cp $output/vcf/$basename.vcf ../bin
        Rscript ../bin/Genome_variant_vis.R
        rm ../bin/*vcf
        mv *.pdf $output/figures

    fi
    echo "running snpEff.."
    [ ! -f $output/vcf/$basename.vcf ] && cp $filepath $output/vcf/$basename.vcf
    snpEff ann NC_045512.2 $output/vcf/$basename.vcf 1> $output/vcf/${basename}_snpeff.tabular #SARS-cov2 isolate Wuhan-Hu-1, complete genome as ref
    [ -f $output/vcf/$basename.vcf ] && rm $output/vcf/$basename.vcf

    echo "running Heatmap.R .."
    cp $output/vcf/${basename}_snpeff.tabular ../bin
    Rscript ../bin/Heatmap.R
    rm ../bin/*.tabular
    mv *.pdf $output/figures

done
# conda env is deactivated when this sub-shell dies
#conda deactivate