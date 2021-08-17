vlnPlotAgeCol <- function(feature,
                          seurat,
                          noClustLab = T,
                          italic = T, 
                          legend = F,
                          star = NULL,
                          markerLab=NULL,
                          featureName = "gene",
                          group = "numclust",
                          ns = "",
                          gap = 0.02,
                          gapSig = 0.02,
                          coordFlip = T,
                          colorAge) {
  
  
  seurat@meta.data$AGE <- factor(seurat@meta.data$AGE,levels = c("Young","Aged"))
  plot <- VlnPlot(object = seurat, 
                  features = feature,
                  group.by = group,
                  pt.size = 0,
                  split.by = "AGE",
                  ncol = 1,
                  split.plot = TRUE,
                  cols = colorAge)  + 
    xlab(label = NULL) +
    ylab(label = NULL)  
  
  if(noClustLab) {
    plot <-plot + theme(axis.text.y = element_blank())
  }
  
  if(italic) {
    plot <- plot + labs(title = bquote(~italic(.(feature)))) 
  }
  
  if(!legend) {
    plot <- plot + NoLegend()
  }
  
  if(coordFlip == T) {
    plot <- plot + coord_flip() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) 
  }
  if(!is.null(star)) {
    
    if(!feature %in% colnames(seurat@meta.data)) {
      seurat@meta.data[,c(feature)] <- GetAssayData(seurat,slot = "data")[feature,]
    }
    
    hits <- star[which(star[,featureName] == feature),c(group,"significant")]
    if(dim(hits)[1]!=0) {
      star <- hits
      
      print(dim(star))
      star$ident <- star[,group]
      for (c in levels(seurat@meta.data[,group])) {
        print(c)
        if (!c %in% star$ident) {
          star <- rbind(star,c(c,ns,c))
        }
      }
    } else {
      star <- data.frame(group = levels(seurat@meta.data[,group]),
                         significant = rep("",length( levels(seurat@meta.data[,group]))))
      
      colnames(star) <- c(group,"significant")
      star$ident <- star[,group]
      
    }
    
    
    
    rownames(star) <- star[,group]
    #print(colnames(star)[1])
    
    
    getMaxPerclust <- function(clust,seurat,feature) {
      res <- max(seurat@meta.data[which(seurat@meta.data[,group] == clust),feature])
    }
    
    star$y <- sapply(star[,group],FUN = getMaxPerclust,feature = feature,seurat = seurat)
    if(star$significant[which.max(star$y)] == "*") {
      lim <- c(min(seurat@meta.data[,feature])-gap/2,
               max(seurat@meta.data[,feature])+ gapSig*max(seurat@meta.data[,feature]))
      if(feature %in% c("Junb(+)","Runx3(+)","Rela(+)")) {
        print("adjusting")
        gapSig <- gapSig +0.1
        lim <- c(min(seurat@meta.data[,feature])-gap/2,
                 max(seurat@meta.data[,feature])+ gapSig*max(seurat@meta.data[,feature]))
      }
    } else {
      lim <- c(min(seurat@meta.data[,feature])-gap,
               max(seurat@meta.data[,feature])+gap)
      if(feature == "Mcpt8") {
        lim <- c(min(seurat@meta.data[,feature])-gap,
                 max(seurat@meta.data[,feature])+0.8)
      }
      
    }
    
    
    print(lim)
    plot <- plot + 
      geom_text(inherit.aes = F, 
                data = star,
                aes(x = ident,
                    y = y + 0.1*max(seurat@meta.data[,feature]),
                    label = significant,size = 60)) + 
      ylim(lim)
  }
  
  if(feature == 'CAM') { #For a better visibility we arrange the scale for CAM signature 
    plot <- plot + scale_y_continuous(breaks = c(-0.5,0.5,1.5))
  }
  #  if(feature == 'Svendsen') { #For a better visibility we arrange the scale for Svendsen signature 
  #   plot <- plot + scale_y_continuous(breaks = c(-0.25,0,0.25))
  # }
  if(feature == 'AgingSvendsen') { #For a better visibility we arrange the scale for CAM signature 
    plot <- plot  + scale_y_continuous(breaks = c(-0.4,0,0.5)) + ggtitle("Svendsen")
  }
  
  if(feature %in% c("HLOD","HEM","TMC")) {
    plot <- plot + scale_y_continuous(breaks = c(-0.6,0,0.6))
  }
  
  ## Add color for cluster marked
  
  colorCluster <- rep("black",length(levels(seurat@meta.data[,group])))
  names(colorCluster) <- levels(seurat@meta.data[,group])
  colorLab <- markerLab[markerLab[,featureName]==feature,group]
  colorCluster[names(colorCluster) %in% colorLab] <- "red"
  plot <- plot + theme(axis.text.x = element_text(colour = colorCluster))
  
  return(plot)
  
}


g_legend<-function(a.gplot,direction = "horizontal"){
  a.gplot <- a.gplot +theme(legend.direction=direction)
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]] 
  return(legend)}

