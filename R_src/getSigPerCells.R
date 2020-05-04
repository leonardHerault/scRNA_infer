getSignatures <- function(m,pheno,cellType,data,id_type = "gene_short_name",outdir="./",padj = 0.05) {
  phenoTest <- pheno
  phenoTest[which(phenoTest!=cellType)] <- "othersHSPC"
  cat("Calling differentially expressed genes (DESeq2\n")
  des <- DESeqDataSetFromMatrix(m, as.data.frame(phenoTest), design=formula(~phenoTest))
  dds <- DESeq(des,parallel=T)
  res <- results(dds,contrast=c("phenoTest",cellType,"othersHSPC"))
  resOrdered <- res[order(res$padj),]
  #save.image("image.Rdata")
  out <- as.data.frame(resOrdered)
  #resMLE<- results(dds, addMLE=TRUE)
  res <- results(dds)
  resOrdered <- res[order(res$padj),]
  #save.image("image.Rdata")
  out <- as.data.frame(resOrdered)
  log2.counts <- log2(counts(dds, normalized=TRUE) + 1)
  colnames(log2.counts) <- colnames(m)
  out <- data.frame(out, log2.counts[rownames(out),],data[rownames(out),"Gene.Symbol"])
  colnames(out)[length(colnames(out))] = "gene_short_name"
  out_save <- out[which(out$padj < 0.05 & out$log2FoldChange<0),]
  out_save <- out_save[order(out_save$padj),]
  write.table(out_save,paste(outdir,"/",paste(cellType,"vs",paste(unique(pheno[which(pheno!=cellType)]),collapse = ""),sep = "_"),".tsv",sep =""),sep="\t",quote=F,col.names = NA)
  
  gene_sig <- as.vector(out_save[which(out_save$padj<padj),id_type])
  return(list(genes = gene_sig, table = out))
}



AddSigScoreMonocle <- function(gbm_cds,ntop,sigRes,cellTypeScore) {
  
  seurat <- exportCDS(gbm_cds,"Seurat")

  sigRes <- sigRes[which(sigRes$padj < 0.05),]
  
  print(head(sigRes))

  gene_sig <- rownames(sigRes[which(sigRes$log2FoldChange<0),])[c(1:50)]
  
  print(head(gene_sig))

  seurat <- AddModuleScore(seurat, genes.list = list(gene_sig), genes.pool = NULL, n.bin = 25,
                  seed.use = 1, ctrl.size = 100, use.k = FALSE,
                  random.seed = 1)
  print("calcul ok")
  print(head(seurat@meta.data))
  pData(gbm_cds)$newScore <- seurat@meta.data[,length(colnames(seurat@meta.data))]
  print("ok)")
  colnames(pData(gbm_cds))[length(colnames(pData(gbm_cds)))] <- cellTypeScore
  print("ok2")
  print(plot_cell_trajectory(gbm_cds,color_by=cellTypeScore) + scale_color_gradient( low="grey",high="red"))
  return(gbm_cds)
}



AddSigScore <- function(seurat,ntop,sigRes,scoreName) {
  
  sigRes <- sigRes[which(sigRes$padj < 0.05),]
  
  print(head(sigRes))
  
  gene_sig <- sigRes[which(sigRes$log2FoldChange<0),"gene_short_name"][c(1:50)]
  
  print(head(gene_sig))
  
  seurat <- AddModuleScore(seurat, genes.list = list(gene_sig), genes.pool = NULL, n.bin = 25,
                           seed.use = 1, ctrl.size = 100, use.k = FALSE, enrich.name= scoreName,
                           random.seed = 1)
  print("calcul ok")
  colnames(seurat@meta.data)[length(colnames(seurat@meta.data))] <- scoreName
  print(head(seurat@meta.data))
  return(seurat)
}

getKeggSig <- function(queryIndex) {
  fullName <- queryIndex$NAME
  name <- strsplit(fullName,split = " ")[[1]][1]
  name <- gsub(name,pattern = "-",replacement = "_")
  name <- paste(name,"_Kegg",sep ="")
  print(name)
  result <- data.frame(strsplit(queryIndex$GENE,";"))
  result_t <- t(result)
  result_vector <- as.vector(result_t[,1])
  result_genes <- result_vector[!grepl("^[0-9]{1,}$", result_vector)]
  
  return(list(name = name,genes = result_genes))
}

getKeggSigList <- function(query) {
  results <- lapply(query,getKeggSig)
  names(results) <- unlist(lapply(results,'[[',1))
  return(results)
}


getMicroArraySig <- function(file) {
  arraySigName <- paste(strsplit(file,split="_")[[1]][2],"_Chambers",sep="")
  arraySig <- as.vector(na.omit(read.table(paste(opt$input_microArraySigDir,file,sep="/"),sep =";",header =T,quote = "\"")$Gene.Symbol))
  return(list(name = arraySigName,genes = arraySig))
}

getMicroArraySigList <- function(fileList) {
  results <- lapply(fileList,getMicroArraySig)
  names(results) <- unlist(lapply(results,'[[',1))
  return(results)
}


read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}


getGeneNameCol <- function(sheet) {
  result <- sheet[,"Gene Symbol"]
  result <- result[which(result != "NA")]
  return(result)
}


getMicroArraySigListXls <- function(allsheets) {
  allsheets <- allsheets[-1]
  results <- lapply(allsheets,getGeneNameCol)
  names(results) <-paste(gsub("\\.txt","",names(allsheets)),"_Chambers",sep="")
  return(results)
}





getBM_vector <- function(goTerm,idType = "external_gene_name",mart) {
  result <- getBM(attributes=c("ensembl_gene_id","external_gene_name","name_1006"),
                  filters=c("go"),
                  values=list(goTerm),mart=mart)
  head(result)
  out <- unique(result[,idType])
  name <- gsub(Term(goTerm),pattern=" |-", replacement= "_")
  name <- paste(name,"_GO",sep="")
  print(name)
  return(list(name = name,genes = out))
}

#Function to score cell in seurat with a gene signature

# DEPRECATED Seurat 2
# scoreCells <- function(seurat,signature,outdir,sigName) {
#   #remove old slot
#   if (!is.null(seurat@meta.data[[sigName]])) {
#     seurat@meta.data[[sigName]] <- NULL
#   }
#   seurat <- AddModuleScore(seurat, genes.list = list(signature), genes.pool = NULL, n.bin = 25,
#                            seed.use = 1, ctrl.size = 100, use.k = FALSE, enrich.name= sigName,
#                            random.seed = 1)
#   colnames(seurat@meta.data)[length(seurat@meta.data)] <- sigName
#   #print(head(seurat@meta.data))
#   if(!is.null(outdir)) {
#   #print("plot")
#   png(paste(outdir,"/",sigName,".png",sep =""))
#   FeaturePlot(seurat,features.plot=sigName)
#   dev.off()
#   }
#   return(seurat)
#   
# }

#For seurat3

scoreCells3 <- function(seurat,signature,outdir,sigName) {
  #remove old slot
  if (!is.null(seurat@meta.data[[sigName]])) {
    seurat@meta.data[[sigName]] <- NULL
  }
  seurat <- AddModuleScore(seurat, features = list(signature), pool = NULL, nbin = 25,
                           seed = 1, ctrl = 100, k = FALSE, name= sigName)
  colnames(seurat@meta.data)[length(seurat@meta.data)] <- sigName
  #print(head(seurat@meta.data))
  if(!is.null(outdir)) {
    #print("plot")
    png(paste(outdir,"/",sigName,".png",sep =""))
    plot(FeaturePlot(seurat,features=sigName))
    dev.off()
  }
  return(seurat)
  
}
## Test cell signature

