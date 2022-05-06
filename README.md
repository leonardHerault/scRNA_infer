
# A novel Boolean network inference strategy to model early hematopoiesis aging
***Leonard Herault***

## Abstract
Hematopoietic stem cell (HSC) aging is a multifactoriel event that leads to changes in HSC properties and function. These changes are intrinsically coordinated and affect the early hematopoiesis, involving hematopoietic stem and progenitor cells (HSPCs). The objective of this work is to better understand the mechanisms and factors controlling these changes. We have therefore developed an original strategy to construct a Boolean network of genes explaining the priming and homeostasis of HSCs. Based on our previous scRNA-seq data, we performed an exhaustive analysis of the transcriptional network and identified active transcription modules or regulons along the differentiation trajectory of selected HSPC states. This global view of transcriptional regulation led us to focus on 15 components, 13 selected TFs (Tal1, Fli1, Gata2, Gata1, Zfpm1, Egr1, Junb, Ikzf1, Myc, Cebpa, Bclaf1, Klf1, Spi1) and 2 complexes regulating the ability of HSC to cycle (CDK4/6 - Cyclin D and CIP/KIP). We then defined the connections controlling the differentiation dynamics of HSC states and constructed an influence graph between the TFs involved in the dynamics by mixing observations from our scRNA-seq data and knowledge from the literature. Then, using answer set programming (ASP) and in silico perturbation analysis, we obtained a Boolean model which is the solution of a Boolean satisfiability problem. Finally, perturbation of the model based on age-related changes revealed important regulations, such as the overactivation of Egr1 and Junb or the loss of Cebpa activation by Gata2, which were found to be relevant for the myeloid bias of aged HSC.
Our work shows the efficiency of the combination of manual and systematic methods to elaborate a Boolean model. The developed strategy led to the proposal of new regulatory mechanisms underlying the differentiation bias of aged HSCs, explaining the decreased transcriptional priming of HSCs to all mature cell types except megakaryocytes.




## Availabilty of data

The single-cell RNA-seq data used in our study are available in the Gene Expression Omnibus database under accession code GSE147729. The workflow is based on the results of our previous work [Hérault et al, 2021](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-00955-z) so the workflow from Herault et al, 2021 should be run first. Mouse TF experiment bed files should be downloaded from [citrome database](http://cistrome.org/db/#/) (available under request) and unarchived as input/mouse_factor.

## Analysis and script

**This workflow is still under development until the manuscript submission and some files are deprecated.**

We provide in this repository the snakemake workflow (Snakefile.py) and its configuration (config/scRNA_infer.yml) that we developed to infer a Boolean gene network from our previously analyzed RNA seq data. It uses [pySCENIC](https://pyscenic.readthedocs.io/en/latest/) command line tools to infer regulons activities in our data. We used [Bonesis](https://github.com/bioasp/bonesis.git) to obtain the possible solutions of the model with our dynamic constraints. See the material and method section of our manuscript for more details.
The final jupyter notebook (included in the workflow) with the model analysis in MP semantic with [mpbn](https://github.com/pauleve/mpbn) needs files (provided here):

    *output/Inference/influenceGraph/infGraphTable45.tsv
    *output/Inference/obsDataDis.csv
    *output/Inference/bonesis/possible_final_solutions.p
    
Bonesis environment used in this study can be installed as follow:

    conda env create -f config/bonesis_env2.yml -n bonesis_env
    conda activate bonesis_env
    cd config/
    git clone https://github.com/bioasp/bonesis.git
    cd bonesis
    pip install --user -e .;
    
jupyter notebook can be launched as follow:   

    jupyter nbconvert --ExecutePreprocessor.timeout=1000000 --to HTML --execute report/reportBonesis.ipynb


## Installation

This snakemake workflow works with conda on Linux.
You need first to download and install [conda with python version 3.7](https://docs.conda.io/en/latest/miniconda.html).
Then once you have downloaded the repository, you can create the snakemake environment with:

    conda env create -n snakemake -f config/snakemake_env.yml

### Activate snakemake environment

    conda activate snakemake

### Command line to launch the  workflow 

    snakemake -j 1 \
        -nps Snakefile.py \
        --configfile config/scRNA_infer.yml \
        --use-conda \
        --forceall --rulegraph > dag.dot

### command lines to launch the snakemake workflow on a cluster with slurm
    
    #export LC_ALL=C\
    #export PATH="/shared/home/lherault/bin/miniconda3/bin:$PATH"
    conda activate snakemake
    snakemake -j 7 \
        -ps Snakefile.py \
        --configfile config/scRNA_infer.yml \
        --use-conda \
        --cluster-config config/cluster.yml \
        --cluster "sbatch -A {cluster.account} \
        -p {cluster.partition} \
        -N {cluster.N} \
        -t {cluster.time} \
        --job-name {cluster.name} \
        --mem {cluster.mem} \
        --cpus-per-task {cluster.cpus-per-task}\
        --output {cluster.output} \
        --error {cluster.error}"
