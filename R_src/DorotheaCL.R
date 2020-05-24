# Script to make analyzis of pseudotime with a gbm_cds ordered 
# command line interface

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

# if (!require("dorothea")) {
#   if (!require("dorothea")) {
#     install.packages("devtools",repos = "http://cran.us.r-project.org") 
#   }
#   print('installing dorothea from github')
#   Sys.setenv(TAR = "/bin/tar")
#   devtools::install_github("saezlab/dorothea")
# }

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("dorothea")) {
  BiocManager::install("dorothea")
}

if (!require("dplyr")) {
  install.packages("dplyr",repos = "http://cran.us.r-project.org")
}

if (!require("getopt")) {
  install.packages("getopt",repos = "http://cran.us.r-project.org")
}

if (!require("Seurat")) {
  install.packages("Seurat",repos = "http://cran.us.r-project.org")
}

if (!require("viper")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")

  BiocManager::install("viper")
}



library(dplyr)
library(Seurat)
library(dorothea)
library(viper)
library(getopt)






spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputRDS',  'i', 1, "character", "REQUIRED: seurat object (.rds)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  'cores',     'c',1, "numeric", "Number of cores"
), byrow=TRUE, ncol=5);

opt = getopt(spec)



# if help was asked, print a friendly message
# and exit with a non-zero error code

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$inputRDS)) {
  cat("Gene filtering of a gbm cds ordered. For scenic use, two filters")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

#set default arguments values

if (is.null(opt$outdir)) {
  opt$outdir = "./"
}

dir.create(opt$outdir,recursive = T,showWarnings = F)

seurat <- readRDS(opt$inputRDS)

## We read Dorothea Regulons for Human:
dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))



## We obtain the regulons based on interactions with confidence level A, B and C
# regulon <- dorothea_regulon_mouse %>%
#   dplyr::filter(confidence %in% c("A","B","C")) %>%
#   df2regulon()


  
regulon <- dorothea_regulon_mouse %>%
  dplyr::filter(confidence %in% c("A","B","C"))
## We compute Viper Scores 
seurat <- run_viper(seurat, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = opt$cores, 
                                 verbose = FALSE))

# seurat <- scViper(seurat, regulon, return_assay = TRUE, 
#                 options = list(method = "scale", minsize = 4, eset.filter = FALSE, 
#                                cores = 1, verbose = FALSE),assay_name = DefaultAssay(seurat))

regulonKept <- rownames(GetAssayData(seurat,assay = "dorothea"))
regulonKeptTable <- dorothea_regulon_mouse[which(regulon$tf %in% regulonKept),]

#Write the dataframe for the tf network contruction with only the regulons kept by scViper
write.table(regulonKeptTable,paste(opt$outdir,"/regulonDorotheaKept.tsv",sep = ""),sep = '\t',quote=F)

saveRDS(seurat,paste(opt$outdir,"/seuratDorothea.rds",sep = ''))









