# Recombination
A collection of scripts for the detection and visualisation of recombinants from SARS-CoV-2 sequence data.

# Visualisation
To run the visualisation script, just call the script with a VCF or FASTA file as its first argument from the visual directory, e.g.: 
```
bash visual.sh ../data/SNP_IMSSC2-91-2022-00035.pass.vcf
```

An output directory can be provided as an optional second argument to the script.
In case you call `visual.sh` on a FASTA file the script will try to run `covSonar` on the file to generate a VCF file for visualisation.