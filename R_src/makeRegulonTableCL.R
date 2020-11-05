##------------------------------------------------------------------------------
## L. Herault
##------------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

suppressMessages(library(getopt))
suppressMessages(library(stringr))
suppressMessages(library(plyr))


# Analysis of seurat 3 integration and clusternig workflow

## -----------------------------------------------------------------------------
## Command line args
## -----------------------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "tfTested",         't',1, "character", "TF tested list (eg TF with a motif in the scenic database",
  "regulonJson",       'r',1, "character", "Regulon main Json file",
  "regulonJsonSupp", "s", 1, "character", "Regulon supp Json files (separated by +)",
  "subConditionName", "n", 1, "character", "correspondong names (separated by +)"
  
), byrow=TRUE, ncol=5);



opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$regulonJson)) {
  cat("Perform diff expression between clusters previously obtained with Seurat3, make a summary table of cluster metrics and signature enrichments")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}


if (is.null(opt$logfc_threshold)) {
  opt$logfc_threshold = 0.25
} 
if (is.null(opt$norm_method)) {
  opt$norm_method = "logNorm"
} 

# Printing option before running

for (o in names(opt)) {
  print(paste(o,":", opt[[o]]))
}

# For testing
# setwd("/shared/projects/scRNA_HSPC_Aging/scRNA_infer/")
# 
# opt <- list()
# opt$tfTested <- "output/publicData/mm_mgi_tfs.txt"
# 
# opt$regulonJson <- "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulonsMetaTest.json"
# 
# opt$regulonJsonSupp <- "output/ScenicRNA_multipleRuns_young/cis_target_maskDropouts/aggregatedRegulonsMeta.json+output/ScenicRNA_multipleRuns_old/cis_target_maskDropouts/aggregatedRegulonsMeta.json"
# opt$subConditionName <- "young+old"
# opt$outdir <- 'output/regulonAnalysis/'

dir.create(opt$outdir,showWarnings = F,recursive = T)

jsonSupp <- strsplit(opt$regulonJsonSupp,split = "\\+")[[1]]
subConditionName <- strsplit(opt$subConditionName,split = "\\+")[[1]]

suppData <- cbind(jsonSupp,subConditionName)

###


tf_mouse <- read.table(opt$tfTested)
length(tf_mouse$V1)

tf_tested <- tf_mouse$V1


# Interaction table regulons and their target genes


makeTable <- function(regulonJsonAndName) {
  
  regulonJsonAndName <- unlist(regulonJsonAndName)
  regulonJson <- regulonJsonAndName[1]
  if(!is.na(regulonJsonAndName[2])) {
    suffixCol <- paste0("_",regulonJsonAndName[2])
  } else {
    suffixCol <- ""
  }
  print(suffixCol)
  
  regulons <- RJSONIO::fromJSON(regulonJson)
  
  regulonTable<- matrix(ncol = 6)
  colnames(regulonTable) <- c("regulon","gene","mor",
                              paste0("recoveredTimes",suffixCol),
                              paste0("importanceMean",suffixCol),
                              paste0("importanceSd",suffixCol))
  
  
  for (i in names(regulons)) {
    mor = 1 
    if (endsWith(i,"(-)")) {
      mor <- -1
    }
    for (t in names(regulons[[i]])) {
      importanceMean <- NA
      importanceSd <- NA
      tf <- str_split_fixed(i,"\\(",n=2)[,1]
      if(tf != t) {
        importanceMean <- regulons[[i]][[t]][[2]]
        importanceSd <- regulons[[i]][[t]][[3]]
      }
      add <- unlist(c(tf,t,mor,regulons[[i]][[t]][[1]][1],importanceMean,importanceSd))
      regulonTable<- rbind(regulonTable,add)
    }
  }
  regulonTable <- data.frame(regulonTable[-1,]) ## first line is NA
  regulonTable$interaction <- paste(regulonTable$regulon,regulonTable$gene,regulonTable$mor,sep = "_") 
  rownames(regulonTable) <- regulonTable$interaction
  return(regulonTable)
}


mainRegulonTable <- makeTable(regulonJsonAndName = c(opt$regulonJson, NULL))

if (!is.null(opt$subConditionName)) {
  conditions <- mapply(list, jsonSupp, subConditionName, SIMPLIFY=F)
  suppRegulonTables <- lapply(conditions, makeTable)
  ## write the supp tables
  for (c in (1:length(conditions))) {
  write.table(suppRegulonTables[[c]],paste0(opt$outdir,"/",conditions[[c]][[2]],"RegulonTable.tsv"),sep = '\t')
  }
}

#Add interaction recovered in one or two condition in the main table

for (table in suppRegulonTables) {
  mainRegulonTable <- cbind(mainRegulonTable,table[rownames(mainRegulonTable),c(4:6)])
}

## Write the tables
write.table(mainRegulonTable,paste0(opt$outdir,"/mainRegulonTable.tsv"),sep = '\t')




