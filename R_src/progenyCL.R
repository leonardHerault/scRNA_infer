# Script to make analyzis of pseudotime with a gbm_cds ordered 
# command line interface

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

if (!require("progeny")) devtools::install_github("saezlab/progeny", ref = "master")
if (!require("getopt")) install.packages("getopt",repos = 'http://cran.us.r-project.org')


library(getopt)
library(progeny)



spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputData',  'i', 1, "character", "REQUIRED: seurat normalized expression data (.csv)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)'
), byrow=TRUE, ncol=5);

opt = getopt(spec)



# if help was asked, print a friendly message
# and exit with a non-zero error code

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$inputData)) {
  cat("Gene filtering of a gbm cds ordered. For scenic use, two filters")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

#set default arguments values

if (is.null(opt$outdir)) {
  opt$outdir = "./"
}


data <- read.csv(opt$inputData,row.names = 1)


progeny_scores = progeny(as.matrix(data), scale=TRUE, organism="Mouse", top=500, perm=1)

write.table(progeny_scores,paste(opt$outdir,"/progeny_scores.tsv",sep = ""),quote = F,sep = '\t')







