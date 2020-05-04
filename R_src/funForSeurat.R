firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# DEPRECATED Seurat2
# getClustInfo <- function(clust,signatures,testHyperSig,seurat,groupingName =NULL,markers) {
#   clustInfo <- list()
#   clustInfo$num_cells <- dim(seurat@meta.data[which(seurat@ident==clust),])[1]
#   clustInfo$percent_cells <- clustInfo$num_cells/dim(seurat@meta.data)[1]
#   
#   # if(!is.null(groupingName)) {
#   #   for (sample_name in unique(seurat@meta.data$sampleName)) {
#   #     clustInfo[[paste(sample_name,"percentCells",sep ="")]] <- dim(seurat@meta.data[which(seurat@ident==clust & seurat@meta.data$sampleName==sample_name),])[1]/clustInfo$num_cells
#   #   }
#   # }
#     
#     percentPhases <- table(seurat@meta.data[which(seurat@ident==clust),"phases"])/length(seurat@meta.data[which(seurat@ident==clust),"phases"]) #In fact this is fraction not percentage
#   
#       
#     # clustInfo$percent_G1_G0 <-percentPhases["G1_G0"]
#     # clustInfo$percent_S <- percentPhases["S"]
#     # clustInfo$percent_G2_M <- percentPhases["G2_M"]
#     
#     for (p in c("G1_G0","S","G2_M")) {
#       if (is.element(p,names(percentPhases))) {
#         clustInfo[[p]] <- percentPhases[p] 
#       } else {
#         clustInfo[[p]] <- 0
#       }
#     }
#     
#     # cell types
#     percentCellType <- table(seurat@meta.data[which(seurat@ident==clust),"CellType"])/length(seurat@meta.data[which(seurat@ident==clust),"CellType"])
#     
#     for (t in unique(seurat@meta.data$CellType)) {
#       if (is.element(t,names(percentCellType))) {
#         clustInfo[[t]] <- percentCellType[t] 
#       } else {
#         clustInfo[[t]] <- 0
#       }
#     }
#     
#     if(!is.null(seurat@meta.data$age)) {
#       
#       # cell types
#     percentAge <- table(seurat@meta.data[which(seurat@ident==clust),"age"])/length(seurat@meta.data[which(seurat@ident==clust),"age"])
#     
#     EnrichAge <- (table(seurat@meta.data[which(seurat@ident==clust),"age"]) / ((table(seurat@meta.data[,"age"])/sum(table(seurat@meta.data[,"age"])))*clustInfo$num_cells)) - 1
#     
#     for (a in unique(seurat@meta.data$age)) {
#       if (is.element(a,names(percentAge))) {
#         clustInfo[[a]] <- percentAge[a] 
#         clustInfo[[paste("enrichAge_",a,sep ="")]] <- EnrichAge[a]
#       } else {
#         clustInfo[[a]] <- 0
#       }
#     }
#       
#     }
#     
#     if(!is.null(seurat@meta.data$percentCellTypeAge)) {
#       
#       # cell types by age
#       percentCellTypeAge <- table(seurat@meta.data[which(seurat@ident==clust),"CellTypeAge"])/length(seurat@meta.data[which(seurat@ident==clust),"CellTypeAge"])
#       
#       for (at in unique(seurat@meta.data$CellTypeAge)) {
#         if (is.element(a,names(percentCellTypeAge))) {
#           clustInfo[[at]] <- percentCellTypeAge[at] 
#         } else {
#           clustInfo[[at]] <- 0
#         }
#       }
#       
#     }
#     
#     
#     
#     clustInfo$median_genes_expressed <- median(seurat@meta.data[which(seurat@ident==clust),"numGenesPerCells"])
#     clustInfo$median_nUMI <- median(seurat@meta.data[which(seurat@ident==clust),"Total_mRNAs"])
#     clustInfo$median_percentMitochGenes <- median(seurat@meta.data[which(seurat@ident==clust),"percentMito"])
#     
#     clustSig <- lapply(signatures,testHyperSig,seurat,clust,markers = markers) #forgot markers arg
#     
#     clustInfo <- c(clustInfo,clustSig)
#     
#     return(clustInfo)
#     
# }


FindAgingMarkers3 <- function(cluster,hspc.combined) {
  hspc.combined$cluster.AGE <- paste(Idents(object = hspc.combined), hspc.combined$AGE, sep = "_")
  Idents(object = hspc.combined) <- "cluster.AGE"
  agingMarkers <- FindMarkers(object = hspc.combined, 
                              ident.1 = paste(cluster,"_Old",sep = ""), 
                              ident.2 = paste(cluster,"_Young",sep = ""), 
                              verbose = FALSE)
  agingMarkers$Cluster <- paste(cluster,"_Old_up",sep ="")
  agingMarkers$Cluster[which(agingMarkers$avg_logFC < 0)] <- paste(cluster,"_Old_down",sep ="")
  agingMarkers$Gene <- rownames(agingMarkers)
  agingMarkers <- agingMarkers[order(agingMarkers$avg_logFC),]
  agingMarkers <- agingMarkers[which(agingMarkers$p_val_adj < 0.05),]
  return(agingMarkers)
}

# DEPRECATED Seurat2
# FindAgingMarkers <- function(cluster,seurat.combined,outdir = "./",logfc.threshold = 0.25) {
#   age <- unique(seurat.combined@meta.data$age)
#   ident.1 <- paste0(cluster,"_age_",age[1])
#   ident.2 <- paste0(cluster,"_age_",age[2])
#   age_effect <- FindMarkers(seurat.combined, ident.1 = ident.1, ident.2 = ident.2, 
#                             print.bar = FALSE,logfc.threshold = logfc.threshold)
#   dir.create(paste(outdir,"/aging_test_on_cluster_",cluster,sep =""),showWarnings = F)
#   
#   print(cluster)
#   
#   age_effect <- age_effect[which(age_effect$p_val_adj < 0.05),]
#   
#   if(dim(age_effect)[1]!=0) {
#     age_effect <- age_effect[order(age_effect$avg_logFC),]
#     age_effect$Gene <- rownames(age_effect)
#     age_effect$Cluster <- NA
#     age_effect[which(age_effect$avg_logFC > 0),"Cluster"] <- "Old_down"
#     age_effect[which(age_effect$avg_logFC < 0),"Cluster"] <- "Old_up"
#     age_effect$Cluster <- paste(age_effect$Cluster,"_cluster_",cluster,sep="")
#     
#     resDir <- paste(outdir,"/aging_test_on_cluster_",cluster,sep ="")
#     
#     gprofileClustResult <- gProfileAnalysis(deg_clust = age_effect,
#                                             outdir = paste(resDir,"/gProfiler", sep =""),
#                                             background = row.names(seurat.combined@data))
#     
#     # with a specific bg
#     
# 
#     colClustRes <- "numclust"
#     cellInClust <- row.names(seurat.combined@meta.data[which(seurat.combined@meta.data[,colClustRes]==cluster),])
# 
#     
#     subSeurat <- SubsetData(seurat.combined,cells.use = cellInClust)
#     print(is(subSeurat))
#     
#     #get expressed genes in this cluster
#     num.cells <- rowSums(as.matrix(subSeurat@data) > 0)
#     genes.use <- names(num.cells[which(num.cells >= 1)])
#     
#     gprofileClustResult <- gProfileAnalysis(deg_clust = age_effect,
#                                             outdir = paste(resDir,"/gProfilerSpecificBg", sep =""),
#                                             background = genes.use)
#   }
#   
#   #############################################################################################################
#   write.table(age_effect,paste(outdir,"/aging_test_on_cluster_",cluster,"/AgingMarkers_of_clust",cluster,".tsv", sep =""),sep = "\t",quote = F,col.names = NA)
#   return(age_effect)
# }


# DEPRECATED Seurat2
# removeNonExpressedGenes <- function(seurat,minPropCellExp) {
#   pd <- new("AnnotatedDataFrame", data = seurat@meta.data)
#   fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(seurat@data)))
#   rownames(fd) <- fd$gene_short_name
#   
#   gbm_cds <- newCellDataSet(seurat@raw.data,
#                             phenoData = pd,
#                             featureData = fd,
#                             lowerDetectionLimit = 0.1,
#                             expressionFamily = negbinomial.size())
#   
#   gbm_cds <- detectGenes(gbm_cds, min_expr = 0.1)
#   
#   
#   
#   print("remove non expressed genes in the subset (non expressed in at least X% of the cells X user option in monocle dp feature 5% in seurat tutorial 0,1%)")
#   fData(gbm_cds)$use_for_seurat <- fData(gbm_cds)$num_cells_expressed > minPropCellExp * ncol(gbm_cds)
#   
#   gbm_to_seurat <- gbm_cds[fData(gbm_cds)$use_for_seurat==T,]
#   gbm_to_seurat <- gbm_cds
#   
#   # Only needed if ensemble id to lazy to code the test for the moment 
#   #rownames(gbm_to_seurat) <- make.unique(fData(gbm_to_seurat)$gene_short_name,sep = "_") #be careful in diff exp results for gene test enrichment
#   
#   if (is.element("Cluster",colnames(pData(gbm_to_seurat)))) {
#     colnames(pData(gbm_to_seurat))[which(colnames(pData(gbm_to_seurat))=="Cluster")] <- "Cluster_monocle"
#   }
#   ##Convert to seurat
#   seurat <- exportCDS(gbm_to_seurat,"Seurat")
#   
#   return(seurat)
#   
# }


## DEPRECATED Seurat 2 Rename ident 

# renameIdent <- function(seurat, old.ident.name,new.ident.name) {
#   
#   seurat@ident <- plyr::mapvalues(x = seurat@ident, from = old.ident.name, to = new.ident.name)
#   seurat@ident <- factor(x = seurat@ident, levels = new.ident.name)
#   
#   return(seurat)
# }




getAGEPropPerClustBarplot <- function(hspc.combined) {
  clusterAge <- ddply(hspc.combined@meta.data,~numclust + AGE,nrow)
  
  propExpect <- table(hspc.combined@meta.data$AGE)/length(hspc.combined@meta.data$AGE)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$AGE)[1]]]
  
  #clusterAGE$numclust <- factor(v$clusterNature , levels = c(""))
  clusterAge$AGE <- factor(clusterAge$AGE , levels = c("Old","Young")) 
  
  
  AGE <- ggplot(data.frame(clusterAge), aes(fill = AGE,y = V1, x=numclust,levels = "Young","Old")) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(hspc.combined@meta.data$AGE)))))+
    scale_y_continuous(name = "Age (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip() + geom_hline(yintercept = propYoungExp)+
    theme(legend.title=element_blank())
  return(AGE)
  
}


getSamplePropPerClustBarplot <- function(hspc.combined) {
  clustersampleName <- ddply(hspc.combined@meta.data,~numclust + sampleName,nrow)
  
  propExpect <- table(hspc.combined@meta.data$sampleName)/length(hspc.combined@meta.data$sampleName)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$sampleName)[1]]]
  
  #clustersampleName$numclust <- factor(v$clusterNature , levels = c(""))
  #clustersampleName$sampleName <- factor(clustersampleName$predicted , levels = c("")) 
  
  
  sampleName <- ggplot(data.frame(clustersampleName), aes(fill = sampleName,y = V1, x=numclust)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(hspc.combined@meta.data$sampleName)))))+
    scale_y_continuous(name = "Sample (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip()+
    theme(legend.title=element_blank()) 
  return(sampleName)
  
}

getRunDatePropPerClustBarplot <- function(hspc.combined) {
  clusterrunDate <- ddply(hspc.combined@meta.data,~numclust + runDate,nrow)
  
  propExpect <- table(hspc.combined@meta.data$runDate)/length(hspc.combined@meta.data$runDate)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$runDate)[1]]]
  
  #clusterrunDate$numclust <- factor(v$clusterNature , levels = c(""))
  #clusterrunDate$runDate <- factor(clusterrunDate$predicted , levels = c("")) 
  
  
  runDate <- ggplot(data.frame(clusterrunDate), aes(fill = runDate,y = V1, x=numclust))  +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(hspc.combined@meta.data$runDate)))))+
    scale_y_continuous(name = "Run date (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip()+
    theme(legend.title=element_blank())  
  return(runDate)
  
}

getPhasePropPerClustBarplot <- function(hspc.combined) {
  clusterphases <- ddply(hspc.combined@meta.data,~numclust + phases,nrow)
  
  propExpect <- table(hspc.combined@meta.data$phases)/length(hspc.combined@meta.data$phases)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$phases)[1]]]
  
  #clusterphases$numclust <- factor(v$clusterNature , levels = c(""))
  #clusterphases$phases <- factor(clusterphases$predicted , levels = c("")) 
  
  
  phases <- ggplot(data.frame(clusterphases), aes(fill = phases,y = V1, x=numclust)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= rev(hue_pal()(length(unique(hspc.combined@meta.data$phases)))))+
    scale_y_continuous(name = "Cell cycle phase (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip() +
    theme(legend.title=element_blank())  
  return(phases)
  
}


getPredictedPropPerClustBarplot <- function(hspc.combined) {
  clusterpredicted <- ddply(hspc.combined@meta.data,~numclust + predicted,nrow)
  
  propExpect <- table(hspc.combined@meta.data$predicted)/length(hspc.combined@meta.data$predicted)[]
  propYoungExp <- propExpect[[unique(hspc.combined@meta.data$predicted)[1]]]
  
  #clusterpredicted$numclust <- factor(v$clusterNature , levels = c(""))
  #clusterpredicted$predicted <- factor(clusterpredicted$predicted , levels = c("")) 
  
  
  predicted <- ggplot(data.frame(clusterpredicted), aes(fill = predicted,y = V1, x=numclust)) +
    geom_bar( stat="identity", position="fill")+
    scale_fill_manual( values= hue_pal()(length(unique(hspc.combined@meta.data$predicted))))+
    scale_y_continuous(name = "Cell type (%)", labels = c(0,25,50,75,100))+
    ylab(label = "")+xlab(label = "") + coord_flip()+
    theme(legend.title=element_blank())  
  return(predicted)
  
}




getEnrichAge <- function(hspc.combined,clustCol ='clusterName',metaCol = "age") {
  
  table <- table(hspc.combined@meta.data[,metaCol],hspc.combined@meta.data[,clustCol])
  
  #Remove null column in case of reclustering has been made
  
  table <- table[,as.vector(which(colSums(table)>0))]
  
  tablePercent <- prop.table(table,2)
  
  propExpect <- table(hspc.combined@meta.data[,metaCol])/length(hspc.combined@meta.data[,metaCol])
  propExpectAge_1<- propExpect[[unique(hspc.combined@meta.data[,metaCol])[1]]]
  propExpectAge_2<- propExpect[[unique(hspc.combined@meta.data[,metaCol])[2]]]
  phyper <- rep(NA,length(colnames(table)))
  enrich <- rep(NA,length(colnames(table)))
  tablePercent <- rbind(tablePercent,enrich,phyper)
  
  
  for (age in unique(hspc.combined@meta.data[,metaCol])) {
    for (cluster in colnames(table)) {
      if(tablePercent[age,cluster] > propExpect[[age]]) {
        cells_pull_marked <- table[age,as.character(cluster)]
        cells_pull <- as.numeric(colSums(table)[as.character(cluster)])
        cells_marked_all <- rowSums(table)[age]
        all_cells <- length(hspc.combined@meta.data[,metaCol])
        
        
        
        p.value <-  phyper(q=cells_pull_marked -1, 
                           m=cells_marked_all,
                           n=all_cells - cells_marked_all, k= cells_pull, lower.tail=FALSE)
        
        tablePercent["enrich",cluster] <- age
        
        tablePercent["phyper",cluster] <- p.value
        
      }
    }
  }
  return(tablePercent)

}






# Make a summary teble of cluster metrics and signature enrichments

getClustTable <- function(rodriguezSig,markers,signatures,seurat,outdir) {
  clusterNames <- c("C1","C2","C3","Mk","Er","Ba","Neu","Mo1","Mo2", "preDC","preB","preT")
  
  RodriguezClustersSig <- lapply(X= c(1:length(clusterNames)),FUN = read_xlsx,path = rodriguezSig)
  
  names(RodriguezClustersSig) <- clusterNames 
  
  getOnlyPos <- function(clustersSig) {
    clusterSig <- clustersSig[which(clustersSig$log.effect > 0),]
    return(clusterSig)
  }
  
  RodriguezClustersSigPos <- lapply(X= RodriguezClustersSig, getOnlyPos)
  
  
  signaturesRodriguez <- lapply(RodriguezClustersSigPos,"[[",1 )
  
  
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  colnames(markers) <- firstup(colnames(markers))
  
  getClustEnrichForRodriguez <- function(clust,signatures,seurat,markers) {
    clustSig <- lapply(signatures,testHyperSig3,seurat,markers,clust)
    propCellTypesLearned <- table(seurat@meta.data[which(seurat@meta.data$numclust==clust),"predicted"])/length(seurat@meta.data[which(seurat@meta.data$numclust==clust),"predicted"])
    clustInfo <- c(clustSig,propCellTypesLearned)
    return(clustInfo)
  }
  
  clust_list <- lapply(unique(markers$Cluster),getClustEnrichForRodriguez,signature=signaturesRodriguez,seurat =seurat,markers =markers) 
  
  names(clust_list) <- paste("cluster_",unique(markers$Cluster),sep="")
  
  clust_table <- as.data.frame(matrix(unlist(clust_list), nrow=length(unlist(clust_list[1]))))
  colnames(clust_table) <- names(clust_list)
  rownames(clust_table) <- names(clust_list[[1]])
  
  clust_df <- as.data.frame(t(clust_table))
  
  write.csv(clust_df,file = paste(outdir,"/clustInfoRodriguez.csv",sep =""),quote = F)
  
  
  
  #Clusters table
  
  clust_table <- data.frame()
  
  print("begin cluster table")
  
  getClustInfo <- function(clust,signatures,seurat,markers) { 

    clustInfo <- list()
    clustInfo$num_cells <- dim(seurat@meta.data[which(seurat@active.ident==clust),])[1]
    clustInfo$percent_cells <- clustInfo$num_cells/dim(seurat@meta.data)[1]
    percentPhases <- table(seurat@meta.data[which(seurat@active.ident==clust),"phases"])/length(seurat@meta.data[which(seurat@active.ident==clust),"phases"]) #In fact this is fraction not percentage
    
    if(!is.null(seurat@meta.data$predicted)) {
      percentPredicted <- table(seurat@meta.data[which(seurat@active.ident==clust),"predicted"])/length(seurat@meta.data[which(seurat@active.ident==clust),"predicted"]) #In fact this is fraction not percentage
      
      for (p in unique(seurat@meta.data$predicted)) {
        print(p)
        if (is.element(p,names(percentPredicted))) {
          clustInfo[[p]] <- percentPredicted[p] 
        } else {
          clustInfo[[p]] <- 0
        }
      }
      
    }
    
    if(!is.null(seurat@meta.data$AGE)) {
      percentAGE <- table(seurat@meta.data[which(seurat@active.ident==clust),"AGE"])/length(seurat@meta.data[which(seurat@active.ident==clust),"AGE"]) #In fact this is fraction not percentage
      
      for (p in unique(seurat@meta.data$AGE)) {
        print(p)
        if (is.element(p,names(percentAGE))) {
          clustInfo[[p]] <- percentAGE[p] 
        } else {
          clustInfo[[p]] <- 0
        }
      }
      
    }
    
    if(!is.null(seurat@meta.data$sampleName)) {
      percentsampleName <- table(seurat@meta.data[which(seurat@active.ident==clust),"sampleName"])/length(seurat@meta.data[which(seurat@active.ident==clust),"sampleName"]) #In fact this is fraction not percentage
      
      for (p in unique(seurat@meta.data$sampleName)) {
        print(p)
        if (is.element(p,names(percentsampleName))) {
          clustInfo[[p]] <- percentsampleName[p] 
        } else {
          clustInfo[[p]] <- 0
        }
      }
      
    }
    
    
    for (p in c("G1_G0","S","G2_M")) {
      if (is.element(p,names(percentPhases))) {
        clustInfo[[p]] <- percentPhases[p] 
      } else {
        clustInfo[[p]] <- 0
      }
    }
    
    
    clustInfo$median_genes_expressed <- median(seurat@meta.data[which(seurat@active.ident==clust),"numGenesPerCells"])
    clustInfo$median_nUMI <- median(seurat@meta.data[which(seurat@active.ident==clust),"Total_mRNAs"])
    clustInfo$median_percentMitochGenes <- median(seurat@meta.data[which(seurat@active.ident==clust),"percentMito"])
    
    
    clustSig <- lapply(signatures,testHyperSig3,seurat,markers,clust) 
    
    clustInfo <- c(clustInfo,clustSig)
    
    
  }
  
  allSignatures <- c(signatures,signaturesRodriguez)
  
  clust_list <- lapply(levels(unique(seurat@active.ident)),getClustInfo,allSignatures,seurat,markers) 
  
  names(clust_list) <- paste("cluster_",levels(unique(seurat@active.ident)),sep="")
  
  print("clust_list ok")
  saveRDS(clust_list,paste(outdir,"/clust_list_save.rds",sep =""))
  
  clust_table <- as.data.frame(matrix(unlist(clust_list), nrow=length(unlist(clust_list[1]))))
  colnames(clust_table) <- names(clust_list)
  rownames(clust_table) <- names(clust_list[[1]])
  
  print("clust_table ok")
  
  clust_df <- as.data.frame(t(clust_table))
  
  write.table(x = clust_df,file = paste(outdir,"/clusters_table.tsv",sep =""),sep="\t",quote=F,col.names = NA)
  
  return(clust_df)
}

  
