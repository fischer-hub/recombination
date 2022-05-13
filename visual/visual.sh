#!/usr/bin/env bash

### Description:
# Takes a vcf or fasta file as input 1 and output dir as optional input 2.
# Output is a figure showing the gene map of the input.

### Credits:
# Plots are created with R scripts from https://github.com/jonas-fuchs/SARS-CoV-2-analyses .
# VCF files from FASTA are created with covSonar from https://gitlab.com/s.fuchs/covsonar .



input=${1%/} # Remove possible trailing /

if [[ $(basename $(pwd)) != "visual" ]]; then
    echo "wrong dir"
    cd visual/
fi

for filepath in $input/*; do

    filename=${filepath##*/}
    basename=${filename%.*}
    echo "currently processing: ${filepath}"
    # check if outdir was provided and set if possible
    if [[ $2 ]]; then
        output=$2
    else
        output="../output"
    fi

        mkdir -p $output/vcf
        mkdir -p $output/figures

    # check if input is a vcf file
    if [ "${filename##*.}" = "vcf" ]; then
        # script expects vcf files to be in same dir
        cp $filepath ../bin
        echo "running Genome_variant_vis.R .."
        Rscript ../bin/Genome_variant_vis.R $filename > /dev/null 2>&1;
        #rm ../bin/*vcf
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

    fi

done

# get vcf file from covsonar db
python3 covsonar/sonar.py var2vcf --db mydb -o $output/vcf/merge.vcf --cpus 8
gunzip -f $output/vcf/merge.vcf.gz
# set float seperator to .
export LC_NUMERIC="en_US.UTF-8";

# calculate something like the allel frequency from AC/AN bc covsonar doesn't do this
awk '/^#/{print $0}!/^#/{split($8, an1, "AN="); split(an1[2], an2, ";"); split($8,info,"AC="); split(info[2],ac,",");printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t"; printf $8";AF="; for(c in ac){printf c/an2[1] ","}; printf "\t"; for (i=9; i<NF; i++) printf $i "\t"; print $NF}' $output/vcf/merge.vcf | awk '/^#/{print $0}!/^#/{printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t";printf substr($8,1,length($8)-1) "\t";for (i=9; i<NF; i++) printf $i "\t"; print $NF}' > $output/vcf/merge_AF.vcf
mv $output/vcf/merge_AF.vcf $output/vcf/merge.vcf

echo "running Genome_variant_vis.R .."
cp $output/vcf/merge.vcf ../bin
Rscript ../bin/Genome_variant_vis.R merge.vcf #> /dev/null 2>&1;
#rm ../bin/*vcf
mv *.pdf $output/figures

for vcf in ../bin/*.vcf; do
    filename=${vcf##*/}
    basename=${filename%.*}

    if command -v snpEff &> /dev/null; then
        echo "running snpEff.. ${vcf}"
        [ ! -f $output/vcf/${filename} ] && cp $vcf $output/vcf/${filename}
        snpEff ann -classic -ud 0 NC_045512.2 $output/vcf/${filename} 1> $output/vcf/${basename}_snpeff.vcf #SARS-cov2 isolate Wuhan-Hu-1, complete genome as ref, -ud 0 ignore up and downstream effects
        echo "running SnpSift"
        SnpSift extractFields -s "," -e "." $output/vcf/${basename}_snpeff.vcf "CHROM"    "POS" "REF"	"ALT"	"FILTER"	"DP"	"AF"	"EFF[*].EFFECT"	"EFF[*].IMPACT"	"EFF[*].FUNCLASS"   "EFF[*].AA" "EFF[*].GENE" > $output/vcf/${basename}.tabular
        cp $output/vcf/${basename}.tabular ../bin
    else
        echo "couldn't annotate VCF with snpEff because it is not available"
    fi
done

echo "running Heatmap.R .."
#Rscript ../bin/Heatmap.R #> /dev/null 2>&1;
Rscript ../bin/Heatmap.R
#rm ../bin/*.tabular
mv *.pdf $output/figures

rm -f mydb
# conda env is deactivated when this sub-shell dies
# conda deactivate