install_github("dami82/easyPubMed")
install.packages("devtools")
install.packages("qdapRegex")
library(easyPubMed)
library(devtools)
library(qdapRegex)
library(tidyverse)

dami_query <- '"Nat Genet"[Journal] AND 2021[PDAT]'
dami_on_pubmed <- get_pubmed_ids(dami_query)
dami_abstracts_xml <- fetch_pubmed_data(dami_on_pubmed)
dami_abstracts_list <- articles_to_list(dami_abstracts_xml)
data<-NULL
PMID<-NULL
for(i in 1:length(dami_abstracts_list)){
  tmp<-article_to_df(pubmedArticle = dami_abstracts_list[[i]], autofill = TRUE, max_chars = 50)[1,]
  pmid<-unique(custom_grep(dami_abstracts_list[[i]], tag = "PMID", format = "char"))[1]
  data<-rbind(data,tmp)
  PMID<-c(PMID,pmid)
  print(i)
}

input<-data.frame(data,PMID)
head(input)
openxlsx::write.xlsx(input,file="Nature.Genetics.2021.xlsx")
