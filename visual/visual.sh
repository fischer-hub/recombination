#!/usr/bin/env bash

### Description:
# Takes a vcf or fasta file as input 1 and output dir as optional input 2.
# Output is a figure showing the gene map of the input.

### Credits:
# Plots are created with R scripts from https://github.com/jonas-fuchs/SARS-CoV-2-analyses .
# VCF files from FASTA are created with covSonar from https://gitlab.com/s.fuchs/covsonar .

while getopts m:i:o:t: opts; do
   case ${opts} in
      m) MULTI_FASTA=${OPTARG} ;;
      t) METADATA=${OPTARG} ;;
      i) INPUT=${OPTARG} ;;
      o) OUTPUT=${OPTARG} ;;
   esac
done

input=${INPUT%/} # Remove possible trailing /

# check if multi_fasta mode is active
if [ ! -z ${MULTI_FASTA+x} ]; then
    mkdir -p temporary
    # split multi fasta file
    awk -v var="$seq_date" -F "|" '/^>/ {close(F) ; F = "temporary/"substr($1,2)".fasta"} {print >> F}' $MULTI_FASTA
    
    # add date from metadata sheet to filename if provided
    if [ ! -z ${METADATA+x} ]; then
        for file in temporary/*; do
            tmpfilename=${file##*/}
            tmpbasename=${tmpfilename%.*}
            seq_date=$(grep ${tmpbasename} ${METADATA} | awk '{print $6}')
            mv $file temporary/${tmpbasename}_${seq_date}.fasta
        done;
    else
            mv $file temporary/${tmpbasename}.fasta
    fi

    input="temporary"
fi


if [[ $(basename $(pwd)) != "visual" ]]; then
    echo "wrong dir"
    cd visual/
fi

for filepath in $input/*; do

    filename=${filepath##*/}
    basename=${filename%.*}
    echo "currently processing: ${filepath}"
    # check if outdir was provided and set if possible
    if [[ ! -z ${OUTPUT+x} ]]; then
        output=$OUTPUT
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
    # assume correct input is a fasta file
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

        echo "running covsonar (merge).."
        # adding all sequences from input to database 'mydb_merge'
        # using eight cpus (the database file will be created if it does not exist)
        python3 covsonar/sonar.py add -f $filepath --db mydb_merge --cpus 8

        echo "running covsonar (single).."
        # adding all sequences from input to database 'mydb'
        # using eight cpus (the database file will be created if it does not exist)
        python3 covsonar/sonar.py add -f $filepath --db mydb_single --cpus 8

        # get vcf files from covsonar db
        python3 covsonar/sonar.py var2vcf --db mydb_single -o $output/vcf/${basename}.vcf --cpus 8
        gunzip -f $output/vcf/${basename}.vcf.gz
        rm mydb_single

        # set float seperator to .
        export LC_NUMERIC="en_US.UTF-8";

        # calculate something like the allel frequency from AC/AN bc covsonar doesn't do this (should always be 1)
        awk '/^#/{print $0}!/^#/{split($8, an1, "AN="); split(an1[2], an2, ";"); split($8,info,"AC="); split(info[2],ac,",");printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t"; printf $8";AF="; for(c in ac){printf c/an2[1] ","}; printf "\t"; for (i=9; i<NF; i++) printf $i "\t"; print $NF}' $output/vcf/${basename}.vcf | awk '/^#/{print $0}!/^#/{printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t";printf substr($8,1,length($8)-1) "\t";for (i=9; i<NF; i++) printf $i "\t"; print $NF}' > $output/vcf/${basename}_AF.vcf
        mv $output/vcf/${basename}_AF.vcf $output/vcf/${basename}.vcf
    fi

done

# get vcf file from covsonar db
python3 covsonar/sonar.py var2vcf --db mydb_merge -o $output/vcf/merge.vcf --cpus 8
gunzip -f $output/vcf/merge.vcf.gz
# set float seperator to .
export LC_NUMERIC="en_US.UTF-8";

# calculate something like the allel frequency from AC/AN bc covsonar doesn't do this
awk '/^#/{print $0}!/^#/{split($8, an1, "AN="); split(an1[2], an2, ";"); split($8,info,"AC="); split(info[2],ac,",");printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t"; printf $8";AF="; for(c in ac){printf c/an2[1] ","}; printf "\t"; for (i=9; i<NF; i++) printf $i "\t"; print $NF}' $output/vcf/merge.vcf | awk '/^#/{print $0}!/^#/{printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t";printf substr($8,1,length($8)-1) "\t";for (i=9; i<NF; i++) printf $i "\t"; print $NF}' > $output/vcf/merge_AF.vcf
mv $output/vcf/merge_AF.vcf $output/vcf/merge.vcf

#echo "running Genome_variant_vis.R .."
cp $output/vcf/merge.vcf ../bin
cp $output/vcf/*.vcf ../bin
#Rscript ../bin/Genome_variant_vis.R merge.vcf #> /dev/null 2>&1;
#rm ../bin/*vcf
#mv *.pdf $output/figures

for vcf in ../bin/*.vcf; do
    filename=${vcf##*/}
    basename=${filename%.*}

    if command -v snpEff &> /dev/null; then
        echo "running snpEff.. ${vcf}"
        [ ! -f $output/vcf/${filename} ] && cp $vcf $output/vcf/${filename}
        snpEff ann -classic -ud 0 NC_045512.2 $output/vcf/${filename} 1> $output/vcf/${basename}_snpeff.vcf #SARS-cov2 isolate Wuhan-Hu-1, complete genome as ref, -ud 0 ignore up and downstream effects
        echo "running SnpSift.."
        SnpSift extractFields -s "," -e "." $output/vcf/${basename}_snpeff.vcf "CHROM"    "POS" "REF"	"ALT"	"FILTER"	"DP"	"AF"	"EFF[*].EFFECT"	"EFF[*].IMPACT"	"EFF[*].FUNCLASS"   "EFF[*].AA" "EFF[*].GENE" > $output/vcf/${basename}.tabular
        cp $output/vcf/${basename}.tabular ../bin
    else
        echo "couldn't annotate VCF with snpEff because it is not available"
    fi
done

if [ ! -z ${MULTI_FASTA+x} ]; then

    echo "running Heatmap.R (multi) .."
    mv bin/merge.tabular bin/dummy

    if [ ! -z ${METADATA+x} ]; then
        Rscript ../bin/Heatmap.R T /dev/null 2>&1;
    elif
        Rscript ../bin/Heatmap.R F /dev/null 2>&1;
    fi
    rm ../bin/*.tabular
    mv *.pdf $output/figures
fi

echo "running Heatmap.R (merge) .."
mv bin/dummy bin/merge.tabular
Rscript ../bin/Heatmap.R F > /dev/null 2>&1;
rm ../bin/*.tabular
mv *.pdf $output/figures

rm -f mydb_merge
rm -f mydb_single
[ -d temporary ] && rm -rf temporary
rm ../bin/*.vcf
# conda env is deactivated when this sub-shell dies
# conda deactivate