#' @title plotExpressionDESeq2
#' @param dds The DESeq2 object
#' @param geneName The character or vector of gene names
#' @export
plotExpressionDESeq2 <- function(dds,geneName){
  library(DESeq2)
  
  if(! "DESeqTransform" %in% class(dds) ){ ## if rlog transformation has not been run
    tmp <- assay(rlog(dds))
  }else{ ## if rlog transformation has been run
    tmp <- assay(dds)
  }
  
  
  
  
  if(length(geneName) ==1){
    temp <- as.data.frame(t(tmp[geneName,] ) )
  }else{
    temp <- as.data.frame(tmp[geneName,] )
  }
  
  
  colnames(temp) <- paste0(dds$group,1:length(dds$group))
  temp$name <- factor(geneName,levels = geneName)
  temp_melt <- reshape2::melt(temp,id.vars = "name")
  temp_melt$Group <- unique(dds$group)[1]
  for(i in 2:length(dds$group)){
    temp_melt$Group[grep(unique(dds$group)[i],temp_melt$variable)] <- unique(dds$group)[i]
  }
  axis.font <- element_text(face = "bold", color = "black")
  ggplot(temp_melt, aes(x=name,y=value,fill=Group))+geom_boxplot()+labs(x="Gene Symbol",y="Log normalized read counts")+
    theme(axis.title =axis.font, axis.text = axis.font)+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black",size = 1),
                       axis.title.x=element_text(size=20, face="bold", hjust=0.5,family = "arial"),
                       axis.title.y=element_text(size=20, face="bold", vjust=0.5, angle=90,family = "arial"),
                       legend.title=element_text(size = 15,face = "bold"),legend.text = element_text(size = 18, face = "bold",family = "arial"),
                       axis.text.x = element_text(size = 15,face = "bold",family = "arial",colour = "black") ,axis.text.y = element_text(size = 15,face = "bold",family = "arial"),
                       plot.title = element_text(size=22, face="bold", hjust=0.5,vjust=0.5,family = "arial"))+
    ggtitle("Gene expression level")
}