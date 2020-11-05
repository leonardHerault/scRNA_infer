##------------------------------------------------------------------------------
## L. Herault
##------------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

suppressMessages(library(getopt))
suppressMessages(library(stringr))
suppressMessages(library(plyr))
suppressMessages(library(igraph))


# Creation of influence graph from scenic results and network from litterature

## -----------------------------------------------------------------------------
## Command line args
## -----------------------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "tfNode",         't',1, "character", "path to TF.txt to take",
  "regulonTable",    'r', 1, "character", "Path to regulons table",
  "addCellCycle",   'c', 0, "logical", "Add cell cycle genes gathered into CycDCD46, INK4, CIPKIP complexes (default FALSE)",
  "addBiblioNet",    'b', 1, "character", "Path to biblio net table to add to the influence graph (separated by +)",
  "importanceThreshold", 'i', 1, "numeric", "discard any interaction (except self loops with no score) below this importance score (default 0)",
  "recoveredTimeThreshold",  "v", 1, "numeric", "discard any interaction recovered in less than this number of runs (default no filter)",
  "selfLoop",    "s", 0, "logical", "Keep selfLoop from regulon table (default TRUE)"
), byrow=TRUE, ncol=5);



opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$regulonTable)) {
  cat("Create influence graph from regulon table")
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

# set deafault arg
if (is.null(opt$addCellCycle)) {
  opt$addCellCycle <- FALSE
}

if(is.null(opt$selfLoop)) {
  opt$selfLoop <- T
} else {
  opt$selfLoop <- F
}


# Printing option before running

for (o in names(opt)) {
  print(paste(o,":", opt[[o]]))
}

# Loading regulon table
regulonTable <- read.table(opt$regulonTable)

# Loading TF nodes
selectedTF <- as.character(read.table(opt$tfNode)$V1)

# Add cell cycle
if (opt$addCellCycle) {
  
  CDK46CycD <- c("Cdk4","Cdk6","Ccnd1","Ccnd2","Ccnd3")
  CIPKIP <- c("Cdkn1b","Cdkn1a","Cdkn1c")
  INK4 <-c("Cdkn2a","Cdkn2b","Cdkn2c","Cdkn2d")
  
  cellCycleGenes <- c(CDK46CycD,CIPKIP,INK4)
  
  infGraphTable <- regulonTable[which(regulonTable$regulon %in% selectedTF & regulonTable$gene %in% c(as.vector(selectedTF),cellCycleGenes)),]
  
  # group cell cycle genes in the three complexes
  levels(infGraphTable$gene) <- c(levels(infGraphTable$gene),c("CIPKIP","INK4","CDK46CycD"))
  infGraphTable$gene[infGraphTable$gene %in% CIPKIP] <- "CIPKIP" 
  infGraphTable$gene[infGraphTable$gene %in% INK4] <- "INK4" 
  infGraphTable$gene[infGraphTable$gene %in% CDK46CycD] <- "CDK46CycD" 
  
} else {
  
  infGraphTable <- regulonTable[which(regulonTable$regulon %in% selectedTF & regulonTable$gene %in% selectedTF),]
  
}

if (!is.null(opt$importanceThreshold)) {
  infGraphTable <- infGraphTable[which(infGraphTable$importanceMean > opt$importanceThreshold | is.na(infGraphTable$importanceMean)),]
}

if (!is.null(opt$recoveredTimeThreshold)) {
  infGraphTable <- infGraphTable[which(infGraphTable$recoveredTimes >= opt$recoveredTimeThreshold ),]
}

if (!opt$selfLoop) {
  infGraphTable <- infGraphTable[which(infGraphTable$tf != infGraphTable$gene),]
}


## Write only the regulonTable
write.table(infGraphTable,paste0(opt$outdir,"/infGraphRegulonTable.tsv",sep = ""),sep = "\t",row.names = F)

## Plot the regulon net
regulonNet <- graph_from_data_frame(infGraphTable,directed = T)
E(regulonNet)$color <- as.factor(get.edge.attribute(regulonNet,"mor"))
E(regulonNet)$lty <- as.factor(get.edge.attribute(regulonNet,"mor")+2)
E(regulonNet)$width <- as.factor(get.edge.attribute(regulonNet,"mor")+2)

png(paste0(opt$outdir,"/regulonGraph.png"),width = 800,height = 800)
plot(regulonNet)
dev.off()

## Plot the regulon net only young
regulonNetY <- graph_from_data_frame(infGraphTable[!is.na(infGraphTable$recoveredTimes_young),],directed = T)
E(regulonNetY)$color <- as.factor(get.edge.attribute(regulonNetY,"mor"))
E(regulonNetY)$lty <- as.factor(get.edge.attribute(regulonNetY,"mor")+2)
E(regulonNetY)$width <- as.factor(get.edge.attribute(regulonNetY,"mor")+2)

png(paste0(opt$outdir,"/regulonGraphY.png"),width = 800,height = 800)
plot(regulonNetY)
dev.off()

## Plot the regulon net only old
regulonNetO <- graph_from_data_frame(infGraphTable[!is.na(infGraphTable$recoveredTimes_old),],directed = T)
E(regulonNetO)$color <- as.factor(get.edge.attribute(regulonNetO,"mor"))
E(regulonNetO)$lty <- as.factor(get.edge.attribute(regulonNetO,"mor")+2)
E(regulonNetO)$width <- as.factor(get.edge.attribute(regulonNetO,"mor")+2)

png(paste0(opt$outdir,"/regulonGraphO.png"),width = 800,height = 800)
plot(regulonNetO)
dev.off()

if (!is.null(opt$addBiblioNet)) {
  print("Adding biblio net...")
  nets <- strsplit(opt$addBiblioNet,split = "\\+")
  infGraphTable$source <- "ScenicMulti"
  for (n in nets[[1]]) {
    tableNet <- read.table(n)
    print(head(n))
    colnames(tableNet) <- c("tf","mor","target")
    levels(tableNet$mor) <- c(levels(tableNet$mor),-1,1)
    tableNet$mor[which(tableNet$mor == "-|")] <- -1
    tableNet$mor[which(tableNet$mor == "->")] <- 1
    
    for (r in rownames(tableNet)) {
      if(tableNet[r,"mor"] == "-?") {
        #print(bonzanni[r,])
        neg <- tableNet[r,] 
        neg$mor <- -1
        tableNet[r,"mor"] <- 1
        tableNet <- rbind(tableNet,neg)
      }
    }
    
    tableNet$source <- strsplit(n,split="/")[[1]][length(strsplit(n,split="/")[[1]])]
    #tableNet$interaction <- paste(tableNet$tf,tableNet$target,tableNet$mor,sep = "_") 
    colnames(infGraphTable)[colnames(infGraphTable) == "regulon"] <- "tf"
    colnames(infGraphTable)[colnames(infGraphTable) == "gene"] <- "target"
    infGraphTable <- rbind(infGraphTable[,c("tf","target","mor","source")],tableNet)
  }
  
}


## Plot the final net
regulonNet <- graph_from_data_frame(infGraphTable,directed = T)
E(regulonNet)$color <- as.factor(get.edge.attribute(regulonNet,"mor"))
E(regulonNet)$lty <- as.factor(as.numeric(get.edge.attribute(regulonNet,"mor"))+2)
E(regulonNet)$width <- as.factor(as.numeric(get.edge.attribute(regulonNet,"mor"))+2)

png(paste0(opt$outdir,"/infGraph.png"),width = 800,height = 800)
plot(regulonNet,autocurve = T)
dev.off()



## Write the final table
write.table(infGraphTable,paste0(opt$outdir,"/infGraphTable.tsv",sep = ""),sep = "\t",row.names = F)


