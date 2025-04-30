# Diploid genome assembly of human fibroblast cell lines enables clone specific variant calling, improved read mapping, and accurate phasing
This repository contains a Snakemake pipeline to reproduce the manuscript results. After editing the `config/config.yaml` file with the file locations of the sequencing data, run:
```
snakemake --configfile config/config.yaml --cores 1
```
This will generate phased diploid assemblies for the BJ and IMR-90 cell lines. 
