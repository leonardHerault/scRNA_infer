suppressMessages(library(Seurat))
suppressMessages(library(getopt))



spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'inputRDS',  'i', 1, "character", "REQUIRED : seurat data analysed object (.RDS generated by prepare_data.R).",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "slot", "s", 1, "character", "Seurat data slot to store in the matrix (default data)"
), byrow=TRUE, ncol=5);

opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$inputRDS)) {
  cat("Create data matrix (features x cells) from the desired slot of the seurat object")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}


if (is.null(opt$slot)) {
  opt$slot <- "data"
}

seurat <- readRDS(opt$inputRDS)
matrix <- as.matrix(GetAssayData(seurat,slot = "counts"))
write.csv(matrix,file = paste(opt$outdir,"/dataMatrix.csv",sep =""),quote = F)
