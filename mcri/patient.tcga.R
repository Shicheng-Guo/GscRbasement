library("TCGAbiolinks")
# 1) receive pid from TCGAbiolinks
pid<-TCGAbiolinks:::getGDCprojects()$project_id
pid<-pid[grep("TCGA",pid)]
# 1) receive pid from github
pid<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/drug_response/Pid.drugResponse.txt",head=F,sep="\t")
pid<-as.character(pid[,1])

drug2csv<-function(clinical.drug){
  bcr_patient_barcode<-clinical.drug$bcr_patient_barcode
  therapy_types<-clinical.drug$therapy_types
  drug_name<-clinical.drug$drug_name
  measure_of_response<-clinical.drug$measure_of_response
  days_to_drug_therapy_start<-clinical.drug$days_to_drug_therapy_start
  days_to_drug_therapy_end<-clinical.drug$days_to_drug_therapy_end
  therapy_ongoing<-clinical.drug$therapy_ongoing
  new.clinical.drug<-data.frame(bcr_patient_barcode,therapy_types,drug_name,measure_of_response,days_to_drug_therapy_start,days_to_drug_therapy_end,therapy_ongoing)
  return(new.clinical.drug)
}

radiation2csv<-function(clinical.radiation){
  bcr_patient_barcode<-clinical.radiation$bcr_patient_barcode
  radiation_type<-clinical.radiation$radiation_type
  radiation_dosage<-clinical.radiation$radiation_dosage
  radiation_treatment_ongoing<-clinical.radiation$radiation_treatment_ongoing
  measure_of_response<-clinical.radiation$measure_of_response
  days_to_radiation_therapy_start<-clinical.radiation$days_to_radiation_therapy_start
  days_to_radiation_therapy_end<-clinical.radiation$days_to_radiation_therapy_end
  measure_of_response<-clinical.radiation$measure_of_response
  new.clinical.radiation<-data.frame(bcr_patient_barcode,radiation_type,radiation_dosage,radiation_treatment_ongoing,days_to_radiation_therapy_start,days_to_radiation_therapy_end,measure_of_response)
  return(new.clinical.radiation)
}

uniquesam<-function(x){
  data<-x[match(unique(x[,1]),x[,1]),]  
  return(data)
}

rlt<-c()
for(i in pid){
  query <- GDCquery(project=i,data.category = "Clinical",file.type = "xml")
  GDCdownload(query)
  clinical.drug <- GDCprepare_clinic(query,"drug")
  drugResponse<-drug2csv(clinical.drug)
  rlt<-rbind(rlt,drugResponse)
}

radiation<-c()
for(i in 1:length(pid)){
  query <- GDCquery(project =pid[i],data.category = "Clinical", file.type = "xml")
  GDCdownload(query)
  clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
  temp<-unique(radiation2csv(clinical.radiation))
  radiation<-rbind(radiation,temp)
}
write.table(radiation,file="tcga.radiationResponse.txt",sep="\t",col.names = NA,row.names = T,quote=F)


immunotherapy<-subset(rlt,therapy_types=="Immunotherapy")
HormoneTherapy<-subset(rlt,therapy_types=="Hormone Therapy")
TargetedMoleculartherapy<-subset(rlt,therapy_types=="Targeted Molecular therapy")
Ancillary<-subset(rlt,therapy_types=="Ancillary")

dim(immunotherapy)
dim(HormoneTherapy)
dim(TargetedMoleculartherapy)
dim(Ancillary)

write.table(immunotherapy,file="immunotherapy.txt",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(HormoneTherapy,file="HormoneTherapy.txt",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(TargetedMoleculartherapy,file="TargetedMoleculartherapy.txt",sep="\t",col.names=NA,row.names=T,quote=F)
write.table(Ancillary,file="Ancillary.txt",sep="\t",col.names=NA,row.names=T,quote=F)

sam1<-subset(clinical,therapy_types=="Chemotherapy")
sam2<-subset(clinical,therapy_types=="Targeted Molecular therapy")
sam3<-subset(clinical,therapy_types=="Hormone Therapy")
sam4<-subset(clinical,therapy_types=="Immunotherapy")
sam5<-subset(radiation,radiation_therapy=="YES")

sam1<-uniquesam(sam1)
sam2<-uniquesam(sam2)
sam3<-uniquesam(sam3)
sam4<-uniquesam(sam4)
sam5<-uniquesam(sam5)

nrow(sam1)
nrow(sam2)
nrow(sam3)
nrow(sam4)
nrow(sam5)

RVNRatio(sam1)
RVNRatio(sam2)
RVNRatio(sam3)
RVNRatio(sam4)
RVNRatio(sam5)

RVNRatio<-function(x){
  R1<-subset(x,measure_of_response=="Complete Response")
  R2<-subset(x,measure_of_response=="Partial Response")
  R3<-subset(x,measure_of_response=="Stable Disease")
  N1<-subset(x,measure_of_response=="Clinical Progressive Disease")
  N2<-subset(x,measure_of_response=="Radiographic Progressive Disease")
  X<-nrow(R1)+nrow(R2)+nrow(R3)
  Y<-nrow(N1)+nrow(N2)
  out<-c(X,Y,X/Y)
  return(out)
}

write.table(sam1[,1],file="Chemotherapy.txt",sep="\t",col.names = F,row.names = F,quote=F)
write.table(sam2[,1],file="TargetedMoleculartherapy.txt",sep="\t",col.names =F,row.names = T,quote=F)
write.table(sam3[,1],file="HormoneTherapy.txt",sep="\t",col.names = F,row.names = F,quote=F)
write.table(sam4[,1],file="Immunotherapy.txt",sep="\t",col.names = F,row.names = F,quote=F)
write.table(sam5[,1],file="radiationtherapy.txt",sep="\t",col.names = F,row.names = F,quote=F)


  
