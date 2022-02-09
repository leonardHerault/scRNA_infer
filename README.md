
# Single-cell RNA-seq assisted synthesis of a Boolean network to model early hematopoiesis aging
***Leonard Herault***

## Abstract
We previously analyzed 15 000 transcriptomes of mouse hematopoietic stem and
progenitor cells (HSPCs) from young and aged mice and characterized the early
differentiation of the hematopoietic stem cells (HSCs) according to age, thanks to cell
clustering and pseudotime analysis 1

. In this study, we propose an original strategy to build
a Boolean gene network explaining HSC priming and homeostasis based on our previous
single cell data analysis and the actual knowledge of these biological processes (graphical
abstract).
We first made an exhaustive analysis of the transcriptional network on selected HSPC
states in the differentiation trajectory of HSCs by identifying regulons, modules formed by a
transcription factor (TFs) and its targets, from the scRNA-seq data., From this global view
of transcriptional regulation in early hematopoiesis, we chose to focus on 15 components,
13 selected TFs (Tal1, Fli1, Gata2, Gata1, Zfpm1, Egr1, Junb, Ikzf1, Myc, Cebpa, Bclaf1,
Klf1, Spi1) and two complexes regulating the ability of HSC to cycle (CDK4/6 - Cyclines D
and CIP/KIP). We then defined the relations in the differentiation dynamics we want to model
((non) reachability, attractors) between the HSPC states that are partial observations of

binarized activity levels of the 15 components. Besides, we defined an influence graph of
possibly involved TF interactions in the dynamic using regulon analysis on our single cell
data and interactions from the literature. Next, using Answer Set Programming (ASP) and
considering these inputs, we obtained a Boolean model as a final solution of a Boolean
satisfiability problem. Finally, we perturbed the model according to aging differences
underlined from our regulon analysis. This led us to propose new regulatory mechanisms at
the origin of the differentiation bias of aged HSCs, explaining the decrease in the
transcriptional priming of HSCs toward all mature cell types except megakaryocytes.



## Availabilty of data

The single-cell RNA-seq data used in our study are available in the Gene Expression Omnibus database under accession code GSE147729.
This workflow rely on the results of our previous work, thus workflow from [Hérault et al, 2021](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-021-00955-z) need to be launched first.
Mouse TF experiment bed files need to be download from [citrome database](http://cistrome.org/db/#/) (available under request) and untar as input/mouse_factor.

## Analysis and script

**This workflow is still under development until the manuscript submission and some files are deprecated.**

We provide in this repository the snakemake workflow (Snakefile.py) and its configuration (config/scRNA_infer.yml) we developped to infer a gene boolean network from our previous analysed the single cell RNA seq data.
It uses [pySCENIC](https://pyscenic.readthedocs.io/en/latest/) command line tools to infer regulons activities in our data. We used [Bonesis](https://github.com/bioasp/bonesis.git) to obtain the possible with out dynamical constraints 

See the material and method section of our manuscript for more details.

The final jupyter notebook (included in the workflow) with the model analysis in MP semantic with [mpbn](https://github.com/pauleve/mpbn) need files (provided here):

    *output/Inference/influenceGraph/infGraphTable45.tsv
    *output/Inference/obsDataDis.csv
    *output/Inference/bonesis/possible_final_solutions.p
    
It can be launched as follow:   

    jupyter nbconvert --ExecutePreprocessor.timeout=1000000 --to HTML --execute report/reportBonesis.ipynb


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
        --use-conda \
        --forceall --rulegraph > dag.dot

### command line to launch the snakemake workflow on a cluster with slurm
    
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
