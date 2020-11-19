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

