#Functions to compute difference of feature in seurat object (because Seurat compute logFC)
#Functions to compute difference of feature in seurat object (because Seurat compute logFC)

getExtAvgDiff <- function(row) {
  res <- NA
  if(row[1] > 0 & row[2] > 0) {
    res <- min(c(row))
  }
  if(row[1] < 0 & row[2] < 0) {
    res <- max(c(row))
  }
  return(res)
}

featureDiff <- function(seurat,cells.1,cells.2,feature) {
  data <- GetAssayData(seurat,slot = "data")
  total.diff <- mean(data[feature,cells.1]) - mean(data[feature,cells.2])
  return(total.diff)
} 


getTrueDiffAging <- function (seurat,table,colIdent = "State",suffix = "") {
  
  table[,paste("avg_diff_",suffix,sep = "")] <- NA
  
  for (rm in rownames(table)) {
    feature <- table[rm,"Gene"]
    cells.1 <- colnames(seurat)[which(seurat[[colIdent]] == table[rm,colIdent] & seurat$AGE == "Old")]
    cells.2 <- colnames(seurat)[which(seurat[[colIdent]] == table[rm,colIdent] & seurat$AGE == "Young")]
    table[rm,paste("avg_diff_",suffix,sep = "")] <- featureDiff(seurat,cells.1,cells.2,feature)
  }
  
  return(table)
}

getTrueDiff <- function (seurat,table,colIdent = "State",suffix = "",colCond = "AGE") {
  if (!is.null(levels(seurat@meta.data[,colCond]))) {
    cond1 = levels(seurat@meta.data[,colCond])[1]
    cond2 = levels(seurat@meta.data[,colCond])[2]
  } else {
    cond1 = unique(seurat@meta.data[,colCond])[1]
    cond2 = unique(seurat@meta.data[,colCond])[2]
  }
  
  table[,paste("avg_diff_",suffix,sep = "")] <- NA
  
  for (rm in rownames(table)) {
    feature <- table[rm,"Gene"]
    cells.1 <- colnames(seurat)[which(seurat[[colIdent]] == table[rm,colIdent] & seurat[[colCond]] == cond1)]
    cells.2 <- colnames(seurat)[which(seurat[[colIdent]] == table[rm,colIdent] & seurat[[colCond]] == cond2)]
    table[rm,paste("avg_diff_",suffix,sep = "")] <- featureDiff(seurat,cells.1,cells.2,feature)
  }
  
  return(table)
  
}




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
  Idents(object = hspc.combined) <- "cluster.AGE"
  print(paste0(cluster,"_Aged"))
  paste0(cluster,"_Young")
  
  agingMarkers <- FindConservedMarkers(hspc.combined,
                                       assay = DefaultAssay(hspc.combined),
                                       grouping.var = "platform",
                                       test.use=test.use,
                                       ident.1 = paste0(cluster,"_Aged"),
                                       ident.2 = paste0(cluster,"_Young"),
                                       pseudocount.use = pseudocount.use,
                                       min.pct = min.pct,
                                       logfc.threshold = logfc.threshold,
                                       max.cells.per.ident = max.cells.per.ident
  )
  
  if(packageVersion("Seurat")> 4) {
    print("renaming columns as in seurat3")
    colnames(agingMarkers)[endsWith(x = colnames(agingMarkers),suffix = "FC")] <- gsub(x = colnames(agingMarkers)[endsWith(x = colnames(agingMarkers),suffix = "FC")],pattern = "log2",replacement = "log" )
  }
  
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
  
  agingMarkers$Cluster <- paste(cluster,"_Aged_up",sep ="")
  agingMarkers$Cluster[which(agingMarkers$min_avg_logFC < 0)] <- paste(cluster,"_Aged_down",sep ="")
  agingMarkers$Gene <- rownames(agingMarkers)
  agingMarkers <- agingMarkers[order(agingMarkers$min_avg_logFC),]
  if (filterOnPadj) {
    agingMarkers <- agingMarkers[which(agingMarkers$A_p_val_adj < 0.05 & agingMarkers$B_p_val_adj < 0.05),]
  }
  return(agingMarkers)
  
}

