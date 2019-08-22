data<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/pheWAS/master/UK_biobank_top_hits_2019.tsv",head=T,sep="\t")
data<-read.table("/home/guosa/hpc/pheWAS/UK_biobank_top_hits.tsv",head=T,sep="\t")
epifactor<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/CIBM/epivar/epifactors.txt")

head(epifactor)
head(data)

newdata<-data[na.omit(match(epifactor[,1],data$nearest_genes)),]
head(newdata)
sort(table(newdata$phenostring))
sort(table(newdata$nearest_genes))

newdata<-data[grep("DNMT",data$nearest_genes),]
write.table(newdata,file="DNMT.pheWAS.txt",sep="\t",quote=F,col.names = NA,row.names = T)
