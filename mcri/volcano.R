res <- read.table("https://raw.githubusercontent.com/Shicheng-Guo/PANC/master/results.txt", header=TRUE)
head(res)
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-1,1)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>0.5), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>0.75), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)>0.75), textxy(log2FoldChange, -log10(pvalue), labs=name, cex=.8))
