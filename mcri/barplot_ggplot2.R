library("ggplot2")
library("extrafont")
input<-data.frame(item=c("DMSO","TNF","LLDT8"),mean=c(14.7,20.51,31.93),sd=c(4.02,4.25,4.25))
input$item <- factor(input$item,levels = c("DMSO", "TNF", "LLDT8"))
p<- ggplot(input, aes(x=item, y=mean, fill=item)) 
p<- p+ geom_bar(stat="identity", color="black",position=position_dodge(),size=1.1) 
p<- p+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9),size=1.1) 
p<- p+labs(title="", x="", y = "Apoptosis Ratio (TUNEL)")
p<- p+ theme_classic()
p<-p+theme(text = element_text(size=25,family="TT Arial"))
p<- p+scale_fill_manual(values=c('red','green','blue'))
print(p)


input<-data.frame(CON=c("0nm","6nm","13nm","25nm","50nm","100nm"),mean=c(100,95.61,92.26,90.56,87.94,86.53),sd=c(0,26.67,20.78,20.85,25.33,27.52))
input$CON <- factor(input$CON,levels = c("0nm","6nm","13nm","25nm","50nm","100nm"))
p<- ggplot(input, aes(x=CON, y=mean, fill=CON)) 
p<- p+ geom_bar(stat="identity", color="black",position=position_dodge(),size=1.1) 
p<- p+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9),size=1.1) 
p<- p+labs(title="", x="", y = "Relative Absorbance at 450nm")
p<- p+ theme_classic()
p<-p+theme(text = element_text(size=25,family="TT Arial"))
p<- p+scale_fill_manual(values=2:7)
print(p)


input<-data.frame(CON=c("0nm","20nm","50nm","100nm","200nm","500nm","1000nm"),mean=c(100,95.61,92.26,90.56,87.94,86.53),sd=c(0,26.67,20.78,20.85,25.33,27.52))
input$CON <- factor(input$CON,levels = c("0nm","6nm","13nm","25nm","50nm","100nm"))
p<- ggplot(input, aes(x=CON, y=mean, fill=CON)) 
p<- p+ geom_bar(stat="identity", color="black",position=position_dodge(),size=1.1) 
p<- p+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9),size=1.1) 
p<- p+labs(title="", x="", y = "Relative Absorbance at 450nm")
p<- p+ theme_classic()
p<-p+theme(text = element_text(size=25,family="TT Arial"))
p<- p+scale_fill_manual(values=2:7)
print(p)


library("GEOquery")
GSE45665 <- getGEO("GSE45665")
data1 <- as.data.frame(exprs(GSE45665[[1]]))

library("GEOquery")
GSE65908 <- getGEO("GSE65908")
data2 <- as.data.frame(exprs(GSE65908[[1]]))

library("GEOquery")
GSE84074 <- getGEO("GSE84074")
data3 <- as.data.frame(exprs(GSE84074[[1]]))

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/RA/lncRNA")
data<-read.table("GSE84074.txt",head=T,sep="\t",row.names = 1)
GPL<-read.table("GPL19640.txt",head=T,sep="\t",row.names = 1)

Ref<-data.frame(rownames(data),GPL[match(rownames(data),rownames(GPL)),])
write.table(Ref,file="GPL19640.Annotation.txt",quote=F,col.names = NA,row.names = T)

gene<-c("MMP1","MMP3","CCL3","CCL4","CCL26","IL6","IL8","IL19","IL24") 
P<-c()
OP<-c()
for(j in match(gene,Ref$TargetID)){
for(i in 1:100){
inr<-sample(c(1,3,5,7,9),10,replace = T)
if(length(unique(inr))>2){
inj<-inr+1
p<-t.test(data[j,inr],data[j,inj])$p.value
P<-c(P,p)
}
}
OP<-c(OP,sum(P<0.05)/length(P))
}
names(OP)<-gene
par(cex.lab=1.5,cex.axis=1.5)
barplot(OP,col=rainbow(9)[1:9],ylab="Power",main="Sample Size=10")




