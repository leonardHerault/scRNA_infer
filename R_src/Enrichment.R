# function to make gene enrichment analyzis with gProfileR

## function to perform gprofler enrcihment analysis on each gene cluster of 
## a differentially expressed gene result (eg data frame of DE genes with a column cluster)
getGeneClustGprofile <- function(deg_clust,
                                 background,
                                 organism = "mmusculus",
                                 hier_filtering = "none",
                                 ordered_query = F) {
  
  resEnrichTest <- lapply(split(as.vector(deg_clust$Gene),f= deg_clust$Cluster),gprofiler,
                    organism = "mmusculus",
                    custom_bg = background,
                    ordered_query = ordered_query,
                    hier_filtering = hier_filtering)
  return(resEnrichTest)
}


## function to plot gprofiler results
plotWriteResEnrich <- function(resEnrichTest,
                          sources=c("BP","keg"),
                          outdir = "./gprofiler",
                          clusterLabel =names(resEnrichTest),
                          colors = brewer.pal(length(resEnrichTest)+1, "Set1")) {
  
  dir.create(outdir,recursive = T,showWarnings = F)
  
  colors <- colors

  for (cluster in clusterLabel) {
    write.table(file= paste(outdir,"/gprofiler_table_clust_",cluster,".tsv", sep = ""),
                            resEnrichTest[[cluster]],quote = F,sep = '\t',row.names = F)
    for (s in sources) {
      results <- resEnrichTest[[cluster]][which(resEnrichTest[[cluster]][,"domain"] == s),]
      results <- results[order(results$p.value),]
      results <- as.data.frame(results[,c("term.id","term.name","p.value")])
      results$logPval <- -log(base = 10,x = results$p.value)
      if(s == "tf") {
        results$term.name <- paste(results$term.name,results$term.id)
      }
      
      ## TO DO : if to much term plot only the first n 
      
      if (length(results$term.name) > 55) {
        png(paste(outdir,"/gprofiler_top30_",cluster,"_",s,".png",sep =""),width = 600,height = 600)
        gp <- ggplot(results[c(1:30),], aes(x=reorder(term.name, logPval), y=logPval)) +
          geom_bar(stat='identity',fill = colors[which(names(resEnrichTest)==cluster)] ) +
          coord_flip() +
          xlab("term name") +
          ylab("-log10(p value)") +
          theme(text = element_text(size = 20))
        print(gp)
        
        dev.off()
        
        png(paste(outdir,"/gprofiler_",cluster,"_",s,".png",sep =""),width = 1400,height = 1400)
        
      } else {
        png(paste(outdir,"/gprofiler_",cluster,"_",s,".png",sep =""),width = 600,height = 600)
      }
      gp <- ggplot(results, aes(x=reorder(term.name, logPval), y=logPval)) +
        geom_bar(stat='identity',fill = colors[which(names(resEnrichTest)==cluster)] ) +
        coord_flip() +
        xlab("term name") +
        ylab("-log10(p value)") +
        theme(text = element_text(size = 20))
      print(gp)
      
      dev.off()
    }
    
  }
}

## function to composate
gProfileAnalysis <- function(deg_clust,
                             background,
                             organism = "mmusculus",
                             hier_filtering = "moderate",
                             ordered_query = F,
                             sources=c("BP","keg","CC",'MF','tf'),
                             outdir = "/gprofiler",
                             clusterLabel =names(resEnrichTest),
                             colors = brewer.pal(length(resEnrichTest)+1, "Set1")) {
  
  resEnrichTest <- getGeneClustGprofile(deg_clust = deg_clust,
                                        background = background,
                                        organism = organism,
                                        hier_filtering = hier_filtering,
                                        ordered_query = ordered_query)
  
  plotWriteResEnrich(resEnrichTest,
                     sources=sources,
                     outdir = outdir,
                     clusterLabel =clusterLabel,
                     colors = colors)
  
  return(resEnrichTest)
  
}

# function for hypergeometric test

# hypergeometric test https://dputhier.github.io/ASG/practicals/go_statistics_td/go_statistics_td_2015.html

#DEPRECATED Seurat 2
# testHyper <- function(gene_set, signature, background,seurat) {
#     signature <- signature[which(is.element(signature,set = rownames(seurat@data)))]
#     genesMarked <- gene_set[which(is.element(gene_set,set = signature))]
#     p.value <-  phyper(q=length(genesMarked) -1, 
#                        m=length(signature),                  
#                        n=length(background) - length(signature), k= length(gene_set), lower.tail=FALSE)
#     return(p.value)
# }

#For Seurat 3
testHyper3 <- function(gene_set, signature, background,seurat) {
  signature <- signature[which(is.element(signature,set = rownames(seurat)))]
  genesMarked <- gene_set[which(is.element(gene_set,set = signature))]
  p.value <-  phyper(q=length(genesMarked) -1, 
                     m=length(signature),                 
                     n=length(background) - length(signature), k= length(gene_set), lower.tail=FALSE)
  return(p.value)
}



testHyperCells <- function(cells_pull, cells_setTarget, allCells) {
  cellsMarked <- gene_set[which(is.element(gene_set,set = cells_set))]
  p.value <-  phyper(q=length(cellsMarked) -1, 
                     m=length(cells_setTarget),
                     n=length(allCells) - length(cells_setTarget), k= length(cells_pull), lower.tail=FALSE)
  return(p.value)
}

## DEPRECATED Seurat2
# testHyperSig <- function(signature,seurat,markers,clust) {  
#   result <- testHyper(gene_set=markers[which(markers$Cluster == clust),"Gene"],
#                       signature=signature,seurat=seurat,background=rownames(seurat@data))
#   return(result)
# }

#For Seurat 3
testHyperSig3 <- function(signature,seurat,markers,clust) { 
  result <- testHyper3(gene_set=markers[which(markers$Cluster == clust),"Gene"],
                       signature=signature,seurat=seurat,background=rownames(seurat))
  return(result)
}




  


