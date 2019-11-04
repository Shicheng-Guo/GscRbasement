args = commandArgs(trailingOnly=TRUE)

snp=as.character(args[1])
locusname=args[2]
chr=args[3]
localhitfile=args[4]
range=as.numeric(args[5])
best.pval=as.numeric(args[6])

# Usage: make.fancy.locus.plot.unix("rs35469986", "TAB1", "22", locus, 10, 0.00005)
# setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/RA/he2020/impute/R3/PDE1A")
# snp="rs7597403"
# locusname="rs7597403"
# chr=2
# localhitfile="rs7597403.local"
# range=6
# best.pval=0.01
# locus <- read.table(localhitfile, header=T, row.names=1)

print(snp)
print(locusname)
print(chr)
print(localhitfile)
print(range)
print(best.pval)

make.fancy.locus.plot <- function(snp, locusname, chr, locus, range, best.pval){
  locus <- read.table(localhitfile, header=T, row.names=1)
  hit <- locus[snp,]
  # size of the region
  #
  min.pos <- min(locus$POS) - 50000
  max.pos <- max(locus$POS) + 50000
  size.pos <- max.pos - min.pos
  center.pos <- min.pos + ( size.pos / 2 )
  center.100kb.pos <- round(center.pos / 100000) * 100000
  offset.100kb.pos <- round((size.pos/3) / 100000) * 100000
  #
  # range of y-axis
  #
  # this dedicates 33% of the yaxis to the genes, labels, recomb rate
  offset <- ( range * 4 / 3 ) - range
  big.range <- range + offset 
  ystart.gene <- - offset
  ystart.recomb <- - offset + (big.range / 8)
  #
  # recombination rate 
  #
  recomb <- read.table(paste("~/hpc/rheumatology/RA/he2020/RecombinationRate/genetic_map_GRCh37_chr", chr, ".txt", sep=""), header=T,check.names = F)
  # recomb <- read.table(paste("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/RA/he2020/RecombinationRate/genetic_map_GRCh37_chr", chr, ".txt", sep=""), header=T,check.names = F)
  keep.recomb <- subset(recomb, recomb[,2] > min.pos & recomb[,2] < max.pos)
  # keep.recomb <- subset(recomb, recomb[,1] > min.pos & recomb[,1] < max.pos)
  #
  # genes in the region
  #
  genelist <- read.table("~/hpc/db/hg19/refGene.hg19.bed", header=F)
  # genelist <- read.table("//mcrfnas2/bigdata/Genetic/Projects/shg047/db/hg19/refGene.hg19.bed", header=F)
  colnames(genelist)<-c("CHR","START","STOP","STRAND","NM","GENE")
  genelist$SIZE=genelist$STOP-genelist$START
  # genelist <- read.table(paste("known_genes_build35_050307_chr", chr, ".txt", sep=""), header=T)
  genes.in.locus <- subset(genelist, (CHR==paste("chr",chr,sep="") & START > min.pos & START < max.pos ) | ( CHR==paste("chr",chr,sep="") & STOP > min.pos & STOP < max.pos) )
  print(genes.in.locus)
  #
  # genotyped markers
  #
  #print(chr)
  markers.in.strong.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.8 & locus$TYPE == "typed"))
  markers.in.moderate.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.5 & locus$RSQR < 0.8 & locus$TYPE == "typed"))
  markers.in.weak.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.2 & locus$RSQR < 0.5 & locus$TYPE == "typed"))
  markers.not.in.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR<0.2 & locus$TYPE == "typed"))
  #markers.typed <- subset(locus, (row.names(locus) != snp & locus$TYPE == "typed"))
  #
  # imputed SNPs
  #
  imputed.in.strong.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.8 & locus$TYPE == "imputed"))
  imputed.in.moderate.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.5 & locus$RSQR < 0.8 & locus$TYPE == "imputed"))
  imputed.in.weak.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.2 & locus$RSQR < 0.5 & locus$TYPE == "imputed"))
  imputed.not.in.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR<0.2 & locus$TYPE == "imputed"))
  #markers.imputed <- subset(locus, (row.names(locus) != snp & locus$TYPE == "imputed"))
  par(mar=c(4,4,3,4))
  #
  # start plot with recombination rate (in background)
  #
  plot(keep.recomb[,2], ystart.recomb + ( ( keep.recomb[,3] / 60 ) * ( 6 * big.range / 8 )), type="l", col="lightblue", lwd=3, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="", main=locusname, axes=F)
  #
  # axes, titles and legends
  #
  mtext(paste("Chromosome", chr, "position (kb)", sep=" "), side=1, line=2.5)
  axis(1, at=c(center.100kb.pos - offset.100kb.pos, center.100kb.pos, center.100kb.pos + offset.100kb.pos), labels=c((center.100kb.pos - offset.100kb.pos) / 1000, center.100kb.pos / 1000, (center.100kb.pos + offset.100kb.pos) / 1000), las=1) 
  axis(2, at=seq(0,range,2), labels=seq(0,range,2), las=1) 
  mtext("Observed (-logP)", side=2, at=(range/2), line=2)
  axis(4, at=c( ystart.recomb, ystart.recomb + (big.range / 4), ystart.recomb + ( 2 * big.range / 4), ystart.recomb + ( 3 * big.range / 4 ) ), labels=c("0","20","40","60"), las=1)
  mtext("Recombination rate (cM/Mb)", side=4, at=(-offset+big.range/2), line=2)
  box()
  lines(c(min.pos, max.pos), c(0,0), lty="dotted", lwd=1, col="black")
  #
  # this is the hit
  #
  points(hit$POS, -(log10(hit$PVAL)), pch=23, cex=2.5, bg="red")
  text(hit$POS, -(log10(hit$PVAL)), labels=c(row.names(hit)), pos=3, offset=1)
  #if ( -(log10(best.pval)) < range ) {
  #  points(hit$POS, -(log10(best.pval)), pch=23, cex=2.5, bg="blue")
  #  text(hit$POS, -(log10(best.pval)), labels=c(paste("P=",best.pval,sep="")), pos=4, offset=2)
  #}else{
  #  points(hit$POS, range, pch=23, cex=2.5, bg="blue")
  #  text(hit$POS, range, labels=c(paste("P=",best.pval,sep="")), pos=4, offset=1)
  #}
  #
  # plot the imputed SNPs
  #
  #points(imputed.not.in.ld$POS, -(log10(imputed.not.in.ld$PVAL)), pch=21, cex=1.0, bg="white")
  #points(imputed.in.weak.ld$POS, -(log10(imputed.in.weak.ld$PVAL)), pch=21, cex=1.0, bg="grey")
  #points(imputed.in.moderate.ld$POS, -(log10(imputed.in.moderate.ld$PVAL)), pch=21, cex=1.5, bg="orange")
  #points(imputed.in.strong.ld$POS, -(log10(imputed.in.strong.ld$PVAL)), pch=21, cex=1.8, bg="red")
  points(imputed.not.in.ld$POS, -(log10(imputed.not.in.ld$PVAL)), pch=23, cex=1.0, bg="grey")
  points(imputed.in.weak.ld$POS, -(log10(imputed.in.weak.ld$PVAL)), pch=23, cex=1.0, bg="grey")
  points(imputed.in.moderate.ld$POS, -(log10(imputed.in.moderate.ld$PVAL)), pch=23, cex=1.0, bg="grey")
  points(imputed.in.strong.ld$POS, -(log10(imputed.in.strong.ld$PVAL)), pch=23, cex=1.0, bg="grey")
  #
  # plot the genotyped markers
  #
  points(markers.not.in.ld$POS, -(log10(markers.not.in.ld$PVAL)), pch=23, cex=1.0, bg="yellow")
  points(markers.in.weak.ld$POS, -(log10(markers.in.weak.ld$PVAL)), pch=23, cex=1.25, bg="orange")
  points(markers.in.moderate.ld$POS, -(log10(markers.in.moderate.ld$PVAL)), pch=23, cex=1.25, bg="purple")
  points(markers.in.strong.ld$POS, -(log10(markers.in.strong.ld$PVAL)), pch=23, cex=1.25, bg="red")
  #points(markers.not.in.ld$POS, -(log10(markers.not.in.ld$PVAL)), pch=23, cex=1.0, bg=hsv(0,markers.not.in.ld$RSQR,1))
  #points(markers.in.weak.ld$POS, -(log10(markers.in.weak.ld$PVAL)), pch=23, cex=1.25, bg=hsv(0,markers.in.weak.ld$RSQR,1))
  #points(markers.in.moderate.ld$POS, -(log10(markers.in.moderate.ld$PVAL)), pch=23, cex=1.25, bg=hsv(0,markers.in.moderate.ld$RSQR,1))
  #points(markers.in.strong.ld$POS, -(log10(markers.in.strong.ld$PVAL)), pch=23, cex=1.25, bg=hsv(0,markers.in.strong.ld$RSQR,1))
  #
  #
  # plot the genes
  #
  for ( i in 1:nrow(genes.in.locus) ) { 
    if ( genes.in.locus[i,]$STRAND == "+" ) {
      arrows(max(genes.in.locus[i,]$START, min.pos), -offset, min(genes.in.locus[i,]$STOP, max.pos), -offset, length=0.05, lwd=2, code=2, lty="solid", col="darkgreen")
    } else {		
      arrows(max(genes.in.locus[i,]$START, min.pos), -offset, min(genes.in.locus[i,]$STOP, max.pos), -offset, length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
    }
    if ( ! is.na(genes.in.locus[i,]$GENE) ) {
      print((genes.in.locus[i,]$STOP+genes.in.locus[i,]$START)/ 2)
      text((genes.in.locus[i,]$STOP+genes.in.locus[i,]$START)/ 2, -offset+( big.range/20 ), labels=genes.in.locus[i,]$GENE, cex=0.8)
  }
  }
}
seed=sample(seq(1,100000,by=1),1)
pdf(paste(snp,"localhit",seed,"pdf",sep="."))
make.fancy.locus.plot(snp, locusname, chr, localhitfile, range, best.pval)
dev.off()



