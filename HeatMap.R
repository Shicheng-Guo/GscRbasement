HeatMap<-function(data,phen,figure="heatmap.pdf",cexRow = 0.01,cexCol = 1.2,Colv=T,Rowv=T){
  library("gplots")
  colnames(data)=phen
  colors <- colorpanel(75,"midnightblue","mediumseagreen","yellow") 
  colors <-bluered(75)
  colors <-greenred(75)
  sidecol<-function(x){
    x<-as.numeric(as.factor(x))
    col<-rainbow(length(table(colnames(data))))
    sapply(x,function(x) col[x])
  }
  ColSideColors=sidecol(phen)
  pdf(figure)
  heatmap.2(data,trace="none",cexRow = cexRow,cexCol = cexCol, ColSideColors=ColSideColors,density.info="none",col=colors,Colv=Colv,Rowv=Rowv,keysize=0.9, margins = c(5, 10))
  dev.off()
}
