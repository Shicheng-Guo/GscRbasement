make.fancy.locus.plot.unix <- function(snp, locusname, chr, locus, range, best.pval) {
  hit <- locus[snp,]
  min.pos <- min(locus$POS) - 10000
  max.pos <- max(locus$POS) + 10000
  size.pos <- max.pos - min.pos
  center.pos <- min.pos + ( size.pos / 2 )
  center.100kb.pos <- round(center.pos / 100000) * 100000
  offset.100kb.pos <- round((size.pos/3) / 100000) * 100000
  offset <- ( range * 4 / 3 ) - range
  big.range <- range + offset 
  ystart.gene <- - offset
  ystart.recomb <- - offset + (big.range / 8)
  recomb <- read.table(paste("~/hpc/rheumatology/RA/he2020/RecombinationRate/genetic_map_GRCh37_chr", chr, ".txt", sep=""), header=T,check.names = F)
  keep.recomb <- subset(recomb, recomb[,2] > min.pos & recomb[,2] < max.pos)
  genelist <- read.table("~/hpc/db/hg19/refGene.hg19.bed", header=F)
  colnames(genelist)<-c("CHR","START","STOP","STRAND","NM","GENE")
  genelist$SIZE=genelist$STOP-genelist$START
  genes.in.locus <- subset(genelist, (CHR==paste("chr",chr,sep="") & genelist$START > min.pos & genelist$START < max.pos ) | ( CHR==paste("chr",chr,sep="") & genelist$STOP > min.pos & genelist$STOP < max.pos) )
  print(genes.in.locus)
  
  markers.in.strong.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.8 & locus$TYPE == "typed"))
  markers.in.moderate.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.5 & locus$RSQR < 0.8 & locus$TYPE == "typed"))
  markers.in.weak.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.2 & locus$RSQR < 0.5 & locus$TYPE == "typed"))
  markers.not.in.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR<0.2 & locus$TYPE == "typed"))
  imputed.in.strong.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.8 & locus$TYPE == "imputed"))
  imputed.in.moderate.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.5 & locus$RSQR < 0.8 & locus$TYPE == "imputed"))
  imputed.in.weak.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR >= 0.2 & locus$RSQR < 0.5 & locus$TYPE == "imputed"))
  imputed.not.in.ld <- subset(locus, (row.names(locus) != snp & locus$RSQR<0.2 & locus$TYPE == "imputed"))
  par(mar=c(4,4,3,4))
  plot(keep.recomb[,2], keep.recomb[,3], type="l", col="lightblue", lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="", main=locusname, axes=F)
  mtext(paste("Chromosome", chr, "position (kb)", sep=" "), side=1, line=2.5)
  axis(1, at=c(center.100kb.pos - offset.100kb.pos, center.100kb.pos, center.100kb.pos + offset.100kb.pos), labels=c((center.100kb.pos - offset.100kb.pos) / 1000, center.100kb.pos / 1000, (center.100kb.pos + offset.100kb.pos) / 1000), las=1) 
  axis(2, at=seq(0,range,2), labels=seq(0,range,2), las=1) 
  mtext("Observed (-logP)", side=2, at=(range/2), line=2)
  axis(4, at=seq(0,big.range,length=4), labels=c("0","20","40","60"), las=1)
  mtext("Recombination rate (cM/Mb)", side=4, at=(-offset+big.range/2), line=2)
  box()
  lines(c(min.pos, max.pos), c(0,0), lty="dotted", lwd=1, col="black")
  points(hit$POS, -(log10(hit$PVAL)), pch=23, cex=2.5, bg="red")
  text(hit$POS, -(log10(hit$PVAL)), labels=c(row.names(hit)), pos=3, offset=1)
  if ( -(log10(best.pval)) < range ) {
    points(hit$POS, -(log10(best.pval)), pch=23, cex=2.5, bg="blue")
    text(hit$POS, -(log10(best.pval)), labels=c(paste("P=",best.pval,sep="")), pos=4, offset=2)
  }else{
    points(hit$POS, range, pch=23, cex=2.5, bg="blue")
    text(hit$POS, range, labels=c(paste("P=",best.pval,sep="")), pos=4, offset=1)
  }
  points(markers.not.in.ld$POS, -(log10(markers.not.in.ld$PVAL)), pch=23, cex=1.0, bg="white")
  points(markers.in.weak.ld$POS, -(log10(markers.in.weak.ld$PVAL)), pch=23, cex=1.25, bg="yellow")
  points(markers.in.moderate.ld$POS, -(log10(markers.in.moderate.ld$PVAL)), pch=23, cex=1.25, bg="orange")
  points(markers.in.strong.ld$POS, -(log10(markers.in.strong.ld$PVAL)), pch=23, cex=1.25, bg="red")
  points(imputed.not.in.ld$POS, -(log10(imputed.not.in.ld$PVAL)), pch=23, cex=1.0, bg="grey")
  points(imputed.in.weak.ld$POS, -(log10(imputed.in.weak.ld$PVAL)), pch=23, cex=1.0, bg="grey")
  points(imputed.in.moderate.ld$POS, -(log10(imputed.in.moderate.ld$PVAL)), pch=23, cex=1.0, bg="grey")
  points(imputed.in.strong.ld$POS, -(log10(imputed.in.strong.ld$PVAL)), pch=23, cex=1.0, bg="grey")
  for ( i in 1:nrow(genes.in.locus)){ 
    if ( genes.in.locus[i,]$STRAND == "+" ) {
      arrows(max(genes.in.locus[i,]$START, min.pos), -offset+i, min(genes.in.locus[i,]$STOP, max.pos), -offset+i, length=0.05, lwd=2, code=2, lty="solid", col="darkgreen")
      text(genes.in.locus[i,]$START + (genes.in.locus[i,]$SIZE/2), -offset+i, labels=genes.in.locus[i,]$GENE, cex=0.8)
      }else{		
      arrows(max(genes.in.locus[i,]$START, min.pos), -offset+i, min(genes.in.locus[i,]$STOP, max.pos), -offset+i, length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
      text(genes.in.locus[i,]$START + (genes.in.locus[i,]$SIZE/2), -offset+i, labels=genes.in.locus[i,]$GENE, cex=0.8)
    }
  }
}

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/RA/he2020/TAB1")
locus <- read.table("local.new", header=T, row.names=1,as.is=T)
make.fancy.locus.plot.unix("rs35469986", "TAB1", "22", locus, 10, 0.00005)
