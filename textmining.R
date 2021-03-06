library(dplyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(readr)
library(tm)

txt<-"//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/HCC/txt"
Reuters <- VCorpus(DirSource(txt, encoding = "UTF-8"),readerControl = list(language = "lat"))
input<-lapply(Reuters[1], as.character)
reuters <- tm_map(Reuters, stripWhitespace)
reuters <- tm_map(reuters, content_transformer(tolower))
reuters <- tm_map(reuters, removeWords, c(stopwords("english")))
dtm <- DocumentTermMatrix(reuters)
findFreqTerms(dtm, 50)
findAssocs(dtm, "methylation", 0.8)
findMostFreqTerms(dtm,n=20)
