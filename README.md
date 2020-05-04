
# Inference of a transcriptionnal factors boolean network driving hematopoietic stem cell fate in young and old HSC

***Leonard Herault***

## Abstract

## Availabilty of data

The single-cell RNA-seq data generated in our study are available in the Gene Expression Omnibus database under accession code GSE147729.
This workflow rely on the results of our previous work.

## Analysis and script

We provide in this repository the snakemake workflow (Snakefile.py) and its configuration (config/scRNA_infer.yml) we developped to infer a transcriptionnal factor boolean network from our previous analysed the single cell RNA seq data.
It uses [pySCENIC](https://pyscenic.readthedocs.io/en/latest/) command line tools and [Dorothea R package][https://github.com/saezlab/dorothea] to infer regulons activities in our data. We also used [progeny][] to infer pathway activities. We used [DCA] to input denoise our single cell data and remove the dropouts.

See the material and method section of our manuscript for more details.
The final step of the workflow produce an html report with the figure shown in our publication.
Our html produced from this [Rmarkdown file](report/final_report.Rmd) can be download [here](report/our_final_report.html) and view in a html browser.


## Installation

This snakemake workflow work with conda on Linux.
You need first to download and install [conda with python version 3.7](https://docs.conda.io/en/latest/miniconda.html).
Then once you have downloaded the repository, you can create the snakemake environment with:

    conda env create -n snakemake -f config/snakemake_env.yml

### Activate snakemake environment

    conda activate snakemake

### Command line to launch the  workflow 

    snakemake -j 1 \
        -nps Snakefile.py \
        --configfile config/scRNA_infer.yml \
        --use-conda
