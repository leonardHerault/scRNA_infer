##------------------------------------------------------------------------------
## L. Herault
##------------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

suppressMessages(library(getopt))
suppressMessages(library(stringr))
suppressMessages(library(plyr))


# Creation of influence graph from scenic results and network from litterature

## -----------------------------------------------------------------------------
## Command line args
## -----------------------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "BetaResults", "b", 1, "character", "Beta result directory"
  
), byrow=TRUE, ncol=5);



opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$BetaResults)) {
  cat("Compilation of Beta results")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}


#For testing
# selectedTF <- c("Gata2","Runx1","Klf1","Cebpa","Gata1","Fli1","Tal1","Spi1","Ikzf1","Junb","Myc","Bclaf1")
# setwd("/shared/projects/scRNA_HSPC_Aging/scRNA_infer/")
# selectedTF <- c("Junb","Stat1","Irf1","Irf9","Myc","Gata2","Spi1","Bclaf1","Cebpa","Tal1","Ikzf1","Fli1","Klf1","Zfpm1","Gata1")
# write.table(selectedTF,"output/ASP/selectedTF.txt",sep = "\t")
# 
# opt <- list()
# opt$outdir <- "output/ASP/"
# opt$regulonTable <- "output/regulonAnalysis_save//mainRegulonTable.tsv"
# opt$addCellCycle <- TRUE
# opt$biblioNet <- "output/Modelling/MatchingDataStateNet.reggraph"
# opt$inportanceThreshold <- NULL
# opt$recoveredTimeThreshold <- 50
# opt$selfLoop <- TRUE

# # set deafault arg
# if (is.null(opt$addCellCycle)) {
#   opt$addCellCycle <- FALSE
# }
# 
# if(is.null(opt$selfLoop)) {
#   opt$selfLoop <- T
# } else {
#   opt$selfLoop <- F
# }


# Printing option before running

for (o in names(opt)) {
  print(paste(o,":", opt[[o]]))
}



## Bone marow cistrome analysis
#mouse_factor <- read.table("../input/mouse_factor_full_QC.txt",header = T,sep = '\t')

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

firstUpThenLow <- function(x) {
  x <- tolower(x)
  x <- firstup(x)
  return(x)
} 

betaFiles <- list.files(opt$BetaResults,pattern = "targets.txt")

cistromeNet <- read.table(paste(opt$BetaResults,betaFiles[1], sep = "/"),sep = "\t",header = F)
colnames(cistromeNet) <- c("chr","start","end","refseqID","score","strand","gene")
regulon <- strsplit(betaFiles[1], "_")[[1]][2]
gsm <- strsplit(betaFiles[1], "_")[[1]][1]

cistromeNet$regulon <- firstUpThenLow(regulon)
cistromeNet$experiment <- gsm

# geneTargets <- matrix(data = NA,
#                       ncol = length(colnames(cistromeNet))+length(colnames(infos)))
# colnames(geneTargets) <- c(colnames(cistromeNet),colnames(infos))

##create dir for graphs

#dir.create(paste(opt$outdir,"/graphs",sep =""), showWarnings=FALSE)


readTargets <- function(f) {
  tmp <- tryCatch(read.table(paste(opt$BetaResults,
                                   f, 
                                   sep = "/"),
                             sep = "\t",
                             header = F), error=function(e) NULL)
}

for (f in betaFiles[-1]) {
  targets <- readTargets(f)
  
  if (!is.null(targets)) {
    regulon <- strsplit(f, "_")[[1]][2]
    gsm <- strsplit(f, "_")[[1]][1]
    
    targets$regulon <- firstUpThenLow(regulon)
    targets$experiment <- gsm
    colnames(targets) <- c("chr","start","end","refseqID","score","strand","gene","regulon","experiment")
    
    
    
    gsm <- strsplit(f, "_")[[1]][1]
    targets$regulon <- firstUpThenLow(strsplit(f, "_")[[1]][2])
    # png(paste(opt$outdir,"/graphs/",strsplit(f,"_t")[[1]][1],".png", sep =""))
    # hist(targets$score,main = strsplit(f,"_t")[[1]][1])
    # dev.off()
    cistromeNet <- rbind(cistromeNet,targets)
  }
}

cistromeTest <- cistromeNet[duplicated(cistromeNet[,c("gene","regulon","experiment")]),]

cistromeNetUnique <- unique(cistromeNet[,c("start","gene","regulon","score","experiment")]) ## be careful one score per experiment and per TSS
cistromeNetUnique$interaction0 <- paste(cistromeNetUnique$regulon,cistromeNetUnique$gene,sep = "_")


## sum TSS scores for each gene to have a score per gene and experiment
cistromeInteraction <- aggregate(cistromeNetUnique$score, by=list(interaction=cistromeNetUnique$interaction0,experiment = cistromeNetUnique$experiment), FUN=sum)
cistromeInteraction[order(cistromeInteraction$interaction),]

cistromeInteractionExpN <- aggregate(cistromeInteraction$experiment,by = list(cistromeInteraction$interaction),FUN = length)

TF_exp <- unique(cistromeNet[,c("regulon","experiment")])
cistromeTFExpN <- table(TF_exp$regulon)

## final cistrome table results with one score per regulation with the number of experiment for it, the number of experiment for the TF
cistromeReg <-  aggregate(cistromeInteraction$x, by=list(interaction=cistromeInteraction$interaction), FUN=sum)

cistromeReg$nbRecovExp <-cistromeInteractionExpN$x
cistromeReg$regulon <- str_split_fixed(cistromeReg$interaction,pattern = "_",n=2)[,1]
cistromeReg$gene <- str_split_fixed(cistromeReg$interaction,pattern = "_",n=2)[,2]
cistromeReg$nbExp <- NA

for (tf in names(cistromeTFExpN)) {
  cistromeReg$nbExp[cistromeReg$regulon == tf] <- cistromeTFExpN[tf]
}

cistromeReg$meanScore <- cistromeReg$x/cistromeReg$nbRecovExp
cistromeReg$adjustedScore <- cistromeReg$x*(cistromeReg$nbRecovExp/cistromeReg$nbExp)
colnames(cistromeReg) <- c("interaction","sumScore","nbRecovExp","regulon","gene","nbExp","meanScore","adjustedScore")

cistromeReg <- cistromeReg[order(cistromeReg$adjustedScore,decreasing = T),c("interaction","regulon","gene","sumScore","nbRecovExp","nbExp","meanScore","adjustedScore")]

write.table(cistromeReg,paste0(opt$outdir,"/cistromeReg.tsv"),quote = F,sep = "\t")




