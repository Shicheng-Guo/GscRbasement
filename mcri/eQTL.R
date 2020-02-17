
file=list.files(pattern="*.egenes.txt")
output<-c()
csel<-c("rs_id_dbSNP151_GRCh38p7","ref","alt","slope","qval")
for(i in 1:length(file)){
  temp<-read.table(file[i],head=T,sep="\t")
  tissuetype=unlist(strsplit(file[i],"[.]"))[1]
  temp<-data.frame(temp[,match(csel,colnames(temp))],SID=tissuetype) 
}