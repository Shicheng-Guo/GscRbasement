BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
res<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/ferroptosis/master/extdata/tcgameta/ferroptosis.tcga.pancancer.chol.overall.rnaseq.dmg.smd.os.hr.csv")
pp<-EnhancedVolcano(res,
                    lab = res$symbol.x,
                    x = 'beta',
                    y = 'pval',
                    xlab="beta",
                    xlim = c(-1, 1),
                    title = 'T vs N',
                    pCutoff = 10e-3,
                    FCcutoff = 0.25,
                    colAlpha = 1,
                    cutoffLineType = 'blank',
                    cutoffLineCol = 'black',
                    cutoffLineWidth = 0.8,
                    hlineCol = c('grey0', 'grey25','grey50','grey75'),
                    hlineType = 'longdash',
                    hlineWidth = 0.8,
                    pointSize=3,
                    labSize=3,
                    gridlines.major = FALSE,
                    gridlines.minor = FALSE)

ggsave(pp,file="ferroptosis.volcano.eps")


library("Haplin")
input<-"ferroptosis"
pvalues=as.numeric(res$pval)
qqplot<-function(pvalues,output=paste(input,"qqplot.eps",sep="")){
  pdf(output)
  pQQ(pvalues, nlabs =length(pvalues), conf = 0.95)
  dev.off()
}
qqplot(pvalues)
