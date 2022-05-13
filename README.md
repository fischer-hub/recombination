# Recombination
A collection of scripts for the detection and visualisation of recombinants from SARS-CoV-2 sequence data.

# Visualisation
To run the visualisation script, just call the script with the directory containing your FASTA or VCF files as its first argument from the visual directory, e.g.: 
```
~/recombination/visual$ ./visual.sh ../data
```

An output directory can be provided as an optional second argument to the script (default dir is `recombination/output`).

In case you call `visual.sh` on a FASTA file(s) the script will try to run `covSonar` on the file(s) to generate a VCF file for visualisation.

Since `covSonar` does not provide the allel frequency field in its VCF output, the allel frequency is estimated from the provided `AC` and `AN` fields from the VCF file. 

The VCF files are then annotated using `snpEff` and formatted with `SnpSift` for later visualization in the `Heatmap.R` script.

NOTE: Both `covSonar` and `snpEff` use the SARS-cov2 isolate Wuhan-Hu-1 complete genome as reference for variant detection.