cminput2bed<-function(cminput,memo){
CHR=paste("chr",cminput$Chromosome,sep="")
START=paste(cminput$Position-1)
END=paste(cminput$Position)
SNP=paste(cminput$SNP)
P=paste(cminput$trait1)
bed<-data.frame(CHR,START,END,SNP,P)
filename=paste(memo,".bed",sep="")
write.table(bed,file=filename,quote=F,col.names=F,row.names=F,sep="\t")
}
