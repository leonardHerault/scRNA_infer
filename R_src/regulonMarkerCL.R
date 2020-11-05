##------------------------------------------------------------------------------
## L. Herault
##------------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## Libraries
## -----------------------------------------------------------------------------

suppressMessages(library(getopt))
suppressMessages(library(stringr))
suppressMessages(library(plyr))
suppressMessages(library(reshape2))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))



source("R_src/computeDiffFun.R")


# Analysis of regulon markers of Seurat cluster

## -----------------------------------------------------------------------------
## Command line args
## -----------------------------------------------------------------------------

spec = matrix(c(
  'help',        'h', 0, "logical",   "Help about the program",
  "inputSeurat",  "i",1, "character", "input seurat object with cluster column names numclust",
  "inputMonocle", "m",1, "character", "input monocle object with state column named State",
  "regulonScore", "r",1, "character", "input regulon score matrix (AUCell)",
  'outdir',     'o',1, "character", 'Outdir path (default ./)',
  "tfNode",         't',1, "character", "path to TF.txt to analyse in details",
  "cellCycle",   'c', 0, "logical", "Add cell cycle genes in the heatmap RNA analysis",
  "scoreDiff", "s",  1, "numeric", "Score differences threshold to select the markers",
  "condition",  "a", 1,  "character", "name of a meta col for a condition diff test per clusters"
), byrow=TRUE, ncol=5);



opt = getopt(spec)


# if help was asked, print a friendly message
# and exit with a non-zero error code

args <- commandArgs()
if ( !is.null(opt$help) | is.null(opt$inputSeurat)) {
  cat("Create influence graph from regulon table")
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}


# For testing
  # setwd("/shared/projects/scRNA_HSPC_Aging/scRNA_infer/")
  # selectedTF <- c("Junb","Stat1","Irf1","Irf9","Myc","Gata2","Gata1","Spi1","Bclaf1","Cebpa","Tal1","Ikzf1","Fli1","Klf1","Zfpm1")
  # write.table(selectedTF,"input/selectedTF.txt",sep = "\t",row.names = F,col.names = F)
# 
# opt <- list()
# opt$outdir <- "output/regulonAnalysis/"
# opt$inputSeurat <- "../herault_et_al/scHSC_herault/report/seurat_report.rds"
# opt$inputMonocle <- "../herault_et_al/scHSC_herault/report/monocle_report.rds"
# opt$tfNode <- "input/selectedTF.txt"
# opt$regulonScore <- "output/ScenicRNA_multipleRuns/AUCell_maskDropouts_save/regulons_enrichment.csv"
# opt$scoreDiff <- 0.002
# opt$cellCycle <- T

# set deafault arg
if (is.null(opt$scoreDiff)) {
  opt$scoreDiff <- 0.002
}

if(is.null(opt$outdir)) {
  opt$outdir <- "./"
}

if(is.null(opt$cellCycle)) {
  opt$cellCycle <- F
}

selectedTF <- as.character(read.table(opt$tfNode)$V1)


# Printing option before running

for (o in names(opt)) {
  print(paste(o,":", opt[[o]]))
}

# loading R objects

seurat <- readRDS(opt$inputSeurat)
#monocle <- readRDS(opt$inputMonocle)

# Loading regulon score table
regulonTable <- read.csv(opt$regulonScore,row.names = 1,check.names = F)

# For the moment only analyse postive regulons

regulonTable <- regulonTable[,colnames(regulonTable)[which(endsWith(colnames(regulonTable),suffix = "+)"))]]

# Adding regulon scores to metadata

# store the data in a new assay in the seurat object of ordered cells
auCell_data <- t(regulonTable)

seurat[["AUCell"]] <- CreateAssayObject(data=auCell_data)


#find cluster regulon markers
Idents(seurat) <- "numclust"

#Set AUCell slot 
DefaultAssay(seurat) <- "AUCell"


clusterRegulon <- FindAllMarkers(seurat,
                               only.pos = T,
                               logfc.threshold= 0,
                               pseudocount.use = 1,
                               min.pct = 0.1)

# Compute true mean difference in score because Seurat comput only logFC
# The featureDiff functions is loaded from ../R_src/computeDiffFun.R file

clusterRegulon$avg_diff <- NA

for (rm in rownames(clusterRegulon)) {
  feature <- clusterRegulon[rm,"gene"]
  cells.1 <- colnames(seurat)[which(seurat$numclust == clusterRegulon[rm,"cluster"])]
  cells.2 <- colnames(seurat)[which(seurat$numclust != clusterRegulon[rm,"cluster"])]
  clusterRegulon[rm,"avg_diff"] <- featureDiff(seurat,cells.1,cells.2,feature)
}

# cut off on p adjusted value
clusterRegulon <- clusterRegulon[which(clusterRegulon$p_val_adj < 0.05 & clusterRegulon$avg_diff > opt$scoreDiff),]


# rename some columns
colnames(clusterRegulon) <- c("p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","regulon","avg_diff")

# order columns
clusterRegulon <- clusterRegulon[,c("p_val","pct.1","pct.2","p_val_adj","cluster","regulon","avg_diff")]

# Average heatmap and binarization per clust


DefaultAssay(seurat) <- "AUCell"

selectedTF <- selectedTF[paste0(c(selectedTF),"(+)")%in% rownames(seurat)]


avgClustAuc <- AverageExpression(
  seurat,
  #slot = "data",
  assays = "AUCell",
  features = paste0(c(selectedTF),"(+)"),
  return.seurat = T
)

png(paste(opt$outdir,"/heatmapAvgAucell.png",sep =""),height = 800,width = 800)
DoHeatmap(avgClustAuc,features = paste0(selectedTF,"(+)")) 
dev.off()

binarizeScaledMatrix <- function(seuratAvg, assay = "RNA") {
  selectedTFMatrix <- GetAssayData(seuratAvg,assay = assay,slot = "scale.data")
  selectedTFMatrixBin <- selectedTFMatrix
  selectedTFMatrixBin[] <- NA
  selectedTFMatrixBin[selectedTFMatrix < -1] <- 0
  selectedTFMatrixBin[selectedTFMatrix > 1] <- 1
  return(selectedTFMatrixBin)
}

selectedTFMatrixBin <- binarizeScaledMatrix(avgClustAuc,assay = "AUCell")


write.table(selectedTFMatrixBin,paste(opt$outdir,"/BinAvgClustAucell.tsv",sep =""),sep = "\t")

# same with RNA level
DefaultAssay(seurat) <- "RNA"

selectedTF <- as.character(read.table(opt$tfNode)$V1)

selectedTF <- selectedTF[selectedTF %in% rownames(seurat)]


CDK46CycDGenes <- c("Cdk4","Cdk6","Ccnd1","Ccnd2","Ccnd3")
CIPKIPGenes  <- c("Cdkn1b","Cdkn1a","Cdkn1c")
INK4Genes  <-c("Cdkn2a","Cdkn2b","Cdkn2c","Cdkn2d")

cellCycleGenes <- c(CDK46CycDGenes ,CIPKIPGenes,INK4Genes)

cellCycleGenes <- cellCycleGenes[cellCycleGenes%in% rownames(seurat)]

if (opt$cellCycle) {
  avgClustRNA <- AverageExpression(
    seurat,
    #slot = "data",
    assays = "RNA",
    features = c(selectedTF,cellCycleGenes),
    return.seurat = T
  )
  
  avgClustRNACC <- data.frame(GetAssayData(avgClustRNA,assay = "RNA",slot = "data"))
  CIPKIP <- colSums(avgClustRNACC[CIPKIPGenes[CIPKIPGenes %in% rownames(avgClustRNACC)],])
  INK4 <- colSums(avgClustRNACC[INK4Genes[INK4Genes %in% rownames(avgClustRNACC)],])
  CDK46CycD <- colSums(avgClustRNACC[CDK46CycDGenes[CDK46CycDGenes %in% rownames(avgClustRNACC)],])
  avgClustRNACC <- rbind(avgClustRNACC,CIPKIP=CIPKIP)
  avgClustRNACC <- rbind(avgClustRNACC,INK4=INK4)
  avgClustRNACC <- rbind(avgClustRNACC,CDK46CycD=CDK46CycD)
  avgClustRNACC <- avgClustRNACC[c(selectedTF,c("CIPKIP","INK4","CDK46CycD")),]
  
  avgClustRNACC_seurat <- CreateSeuratObject(counts = avgClustRNACC)
  avgClustRNACC_seurat <- ScaleData(avgClustRNACC_seurat)
  Idents(avgClustRNACC_seurat) <- colnames(avgClustRNACC_seurat)
  
  
  png(paste(opt$outdir,"/heatmapAvgRNA.png",sep =""),height = 800,width = 800)
  DoHeatmap(avgClustRNACC_seurat,features = c(selectedTF,c("CIPKIP","INK4","CDK46CycD")))
  dev.off()
  
  selectedTFMatrixBinRNA <- binarizeScaledMatrix(avgClustRNACC_seurat,assay = "RNA")
  write.table(selectedTFMatrixBinRNA,paste(opt$outdir,"/BinAvgClustRNA.tsv",sep =""),sep = "\t")
  
  
  
} else {
  png(paste(opt$outdir,"/heatmapAvgRNA.png",sep =""),height = 800,width = 800)
  DoHeatmap(avgClustRNA,features = selectedTF) 
  dev.off()
  
  selectedTFMatrixBinRNA <- binarizeScaledMatrix(avgClustRNA,assay = "RNA")
  write.table(selectedTFMatrixBinRNA,paste(opt$outdir,"/BinAvgClustRNA.tsv",sep =""),sep = "\t")
  
}

rownames(selectedTFMatrixBin) <- selectedTF[paste0(selectedTF,"(+)") %in% rownames(selectedTFMatrixBin)]
selectedTFMatrixBin <- rbind(selectedTFMatrixBin,selectedTFMatrixBinRNA[c(selectedTF[!selectedTF %in% rownames(selectedTFMatrixBin)],"CIPKIP","INK4","CDK46CycD"),])
consensusMatrixBin <- selectedTFMatrixBin
consensusMatrixBin[] <- NA
consensusMatrixBin[which(selectedTFMatrixBin == selectedTFMatrixBinRNA)] <- selectedTFMatrixBin[which(selectedTFMatrixBin == selectedTFMatrixBinRNA)]

consensusMatrixBinHm <- melt(consensusMatrixBin)

write.table(consensusMatrixBin,paste(opt$outdir,"/binMatConsensus.tsv",sep =""),sep = "\t")

png(paste(opt$outdir,"/heatmapConsensus.png",sep =""),height = 800,width = 800)
ggplot(consensusMatrixBinHm,aes(x=Var2, y=Var1, fill=value)) + scale_fill_gradient(low="blue",high = "red") +
  geom_tile()          
dev.off()

###########################################################################################################################

FindAgingMarkers4 <- function(cluster,
                              hspc.combined,
                              max.cells.per.ident = Inf,
                              filterOnPadj = T,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              keepDiverging = F,
                              test.use = "wilcox",
                              identCol = "numclust",
                              pseudocount.use = 1) {
  
  hspc.combined$cluster.AGE <- paste(hspc.combined@meta.data[,identCol],
                                     hspc.combined@meta.data[,"AGE"], sep = "_")
  
  test <- table( hspc.combined$cluster.AGE, hspc.combined$platform)
  
  if( (test[paste0(cluster,"_Old"),"A"]>3&test[paste0(cluster,"_Old"),"B"]>3) &
      (test[paste0(cluster,"_Young"),"A"]>3&test[paste0(cluster,"_Young"),"B"]>3)) {
    
    Idents(object = hspc.combined) <- "cluster.AGE"
    print(paste0(cluster,"_Old"))
    paste0(cluster,"_Young")
    
    
    
    agingMarkers <- FindConservedMarkers(hspc.combined,
                                         assay = DefaultAssay(hspc.combined),
                                         grouping.var = "platform",
                                         test.use=test.use,
                                         ident.1 = paste0(cluster,"_Old"),
                                         ident.2 = paste0(cluster,"_Young"),
                                         pseudocount.use = pseudocount.use,
                                         min.pct = min.pct,
                                         logfc.threshold = logfc.threshold,
                                         max.cells.per.ident = max.cells.per.ident
    )
    print(dim(agingMarkers))
    #Add a column the min abs(logFC) 
    agingMarkers$min_avg_logFC <- NA
    
    if(!keepDiverging) {
      #print("keep only markers with conserved logfc sign")
      agingMarkers <- agingMarkers[which(agingMarkers$B_avg_logFC/agingMarkers$A_avg_logFC > 0),]
    } else {
      #print("combined p_val of diverging markers between two batch are set to 1 and their min_avg_logfc to 0")
      agingMarkers[which(agingMarkers$B_avg_logFC/agingMarkers$A_avg_logFC < 0),"min_avg_logFC"] <- 0
      agingMarkers[which(agingMarkers$B_avg_logFC/agingMarkers$A_avg_logFC < 0),"minimump_p_val"] <- 1
    }
    
    for (g in rownames(agingMarkers)) {
      if(agingMarkers[g,"A_avg_logFC"] < 0) {
        agingMarkers[g,"min_avg_logFC"] <- max(agingMarkers[g,c("A_avg_logFC","B_avg_logFC")])
      } else {
        agingMarkers[g,"min_avg_logFC"] <- min(agingMarkers[g,c("A_avg_logFC","B_avg_logFC")])
      }
    }
    
    agingMarkers$Cluster <- paste(cluster,"_Old_up",sep ="")
    agingMarkers$Cluster[which(agingMarkers$min_avg_logFC < 0)] <- paste(cluster,"_Old_down",sep ="")
    agingMarkers$Gene <- rownames(agingMarkers)
    agingMarkers <- agingMarkers[order(agingMarkers$min_avg_logFC),]
    if (filterOnPadj) {
      agingMarkers <- agingMarkers[which(agingMarkers$A_p_val_adj < 0.05 & agingMarkers$B_p_val_adj < 0.05),]
    }
  }else {
    agingMarkers <- NULL 
  }
  return(agingMarkers)
  
}


if (!is.null(opt$condition)) { ## only implemented for a meta col AGE for the moment
  DefaultAssay(seurat) <- "AUCell"
  
  
  
  agingRegulonPerClust <- lapply(levels(seurat$numclust),
                                 FindAgingMarkers4,
                                 seurat,
                                 identCol = "numclust",
                                 test.use = "wilcox",
                                 keepDiverging= T,
                                 logfc.threshold = 0,
                                 pseudocount.use = 1,
                                 min.pct = 0.1,
                                 filterOnPadj = F)
  
  
  names(agingRegulonPerClust) <- levels(seurat$numclust)
  
  agingRegulonPerClustTable <- agingRegulonPerClust[[1]]
  for (s in levels(seurat$numclust)[-1]) {
    agingRegulonPerClustTable <- rbind(agingRegulonPerClustTable,agingRegulonPerClust[[s]])
  }
  
  dim(agingRegulonPerClustTable)
  
  table(agingRegulonPerClustTable$Cluster)
  
  #Exclude AUCell score differences with a combined pval < 0.05
  agingRegulonPerClustTable <- agingRegulonPerClustTable[which(agingRegulonPerClustTable$minimump_p_val < 0.05),]
  
  dim(agingRegulonPerClustTable)
  
  agingRegulonPerClustTable$numclust <- str_split_fixed(agingRegulonPerClustTable$Cluster,pattern = "_",n=3)[,1]
  
  
  agingRegulonPerClustTable <- getTrueDiffAging(seurat[,which(seurat$platform == "A")],
                                                agingRegulonPerClustTable,colIdent = "numclust",
                                                suffix = "A")
  
  agingRegulonPerClustTable <- getTrueDiffAging(seurat[,which(seurat$platform == "B")],
                                                agingRegulonPerClustTable,colIdent = "numclust",
                                                suffix = "B")
  
  
  agingRegulonPerClustTable$min_avg_diff <- apply(agingRegulonPerClustTable[,c("avg_diff_A","avg_diff_B")], 1, FUN=getExtAvgDiff) 
  
  
  dim(agingRegulonPerClustTable)
  
  #Exclude AUCell score differences with an opposite sign of variation (pval set to NA by FindAgingMarkers)
  agingRegulonPerClustTable <- na.exclude(agingRegulonPerClustTable)
  
  dim(agingRegulonPerClustTable)
  
  
  
  colnames(agingRegulonPerClustTable)[which(colnames(agingRegulonPerClustTable) == "Cluster")] <- "group"
  colnames(agingRegulonPerClustTable)[which(colnames(agingRegulonPerClustTable) == "Gene")] <- "regulon"
  
  
  
  agingRegulonPerClustTable <- agingRegulonPerClustTable[,c("regulon","numclust","group","minimump_p_val","max_pval","min_avg_diff","A_p_val","avg_diff_A","A_pct.1","A_pct.2","B_p_val","avg_diff_B","B_pct.1","B_pct.2")]
  
  
  
  
  # Filter on average differences before writing final table
  agingRegulonPerClustTable <- agingRegulonPerClustTable[which(abs(agingRegulonPerClustTable$min_avg_diff) > opt$scoreDiff),]
  
  write.table(agingRegulonPerClustTable,paste(opt$outdir,"/condDiffRegulonTable.tsv",sep =""),sep = "\t")
}


saveRDS(seurat,paste0(opt$outdir,"/seuratAUC.rds"))

write.table(clusterRegulon,paste(opt$outdir,"/clusterMarkerRegulonTable.tsv",sep =""),sep = "\t")



