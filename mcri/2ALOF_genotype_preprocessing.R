source("/home/yez/Schrodi_2ALOF/rpgm/2ALOF_process_sub.R")
setwd("/mnt/bigdata/Genetic/Projects/EXOMECHIP_MARSHFIELD/tmp")
outpath = "/mnt/bigdata/Genetic/Projects/Schrodi_2ALOF/data/phased_data_preprocess/"
ann_dat = read.csv("FinalRelease_QC_20140311_Team1_Marshfield.10302018.annovar.hg19_multianno.csv",header=TRUE)
sam_dat = read.delim("FinalRelease_QC_Phenotypes_Marshfield_20140224_Team1.txt",header=TRUE)
## define parameters
## chromosomes
chr_vec = c(1:22)
## main program
## preprocess the genotype data to only select the ones are needed
marker_summary = matrix(-1,nrow=length(chr_vec),ncol=2)
for(i in 1:length(chr_vec)){
  ped_in = paste("exomechip_SNV_PASS_BEAGLE_chr",i,"_phased.ped",sep="")
  map_in = paste("exomechip_SNV_PASS_BEAGLE_chr",i,"_phased.map",sep="")
  map_dat = read.table(map_in,header=FALSE)
  snp_num = dim(map_dat)[1] ## total number of markers phased
  pos1=grep("splicing|exonic|UTR",ann_dat$Func.refGene)
  pos2=grep("Name",ann_dat$gwasCatalog)
  pos3=grep("Name",ann_dat$wgRna)
  pos4=grep("Name",ann_dat$targetScanS)
  pos5=grep("Name",ann_dat$tfbsConsSites)
  pos_keep = unique(c(pos1,pos2,pos3,pos4,pos5))
  ann_dat_map4 = ann_dat[pos_keep,]
  vcf_info_dat<-read.table(map_in,head=T,sep="") 
  vcf_info_dat_map2 = vcf_info_dat[pos_keep,]    
  map_dat2 = map_dat[pos_keep,]
  snp_out_info = data.frame(ann_dat_map4,vcf_info_dat_map2,map_dat2)
  marker_out_list = map_dat2$V2
  out_marker = paste(outpath,"chr",i,"_2ALOF_SNP_list",sep="")    
  marker_summary[i,1] = snp_num
  marker_summary[i,2] = length(marker_out_list)
  out_marker2 = paste(outpath,"chr",i,"_2ALOF_SNP_info_comb.RData",sep="")
  save(snp_out_info,file=out_marker2)
  print(i)
}
marker_summary2 = data.frame(CHR=chr_vec,marker_summary)
colnames(marker_summary2)[-1] = c("SNP_TOT","SNP_SEL")
outf = paste(outpath,"2ALOF_marker_summary.csv",sep="")
write.csv(marker_summary2,file=outf,row.names=FALSE)