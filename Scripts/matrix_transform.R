# Title: Data Integration and Transformation
# Author: Aijin Wang, Contact details: aijinwang3@gmail.com
# Script info: Transform raw data into matrices to record gene mutation 
#                       y/n and counts
# 2018 QSURE program at Memorial Sloan Kettering Cancer Center

setwd("~/TCGA/Data/Raw")
library(magrittr)
library(dplyr)
library(tidyr)
library(readxl)

# get somatic gene list ---------------------------------------------------
gene_list = read_excel("gene_list.xlsx")
somatic_gene = gene_list[gene_list$Somatic_drivers == 1,]$HUGO_Symbol


# prepare dataset, match gene list with raw data --------------------------
load("genes.Rdata")
tumor = genes %>% 
  select(Tumor_Sample_Barcode, Hugo_Symbol) %>% 
  filter(Hugo_Symbol %in% somatic_gene) %>%
  arrange(Tumor_Sample_Barcode)


## The above processed data is stored in "tumor.RData" for easy use


# create mutation count table ---------------------------------------------
load("tumor.RData")
somatics_count = tumor %>%
  group_by(Tumor_Sample_Barcode,Hugo_Symbol) %>% 
  tally() %>% 
  spread(Hugo_Symbol,n,fill=0) %>% as.data.frame
rownames(somatics_count) = somatics_count[,1]
somatics_count = somatics_count[,-1]
write.csv(somatics_count, file = "~/TCGA/Data/New/somatic_count.csv")


# create binary mutation count table --------------------------------------
somatics_binary = ifelse(somatics_count>0,1,0)
write.csv(somatics_binary, file = "~/TCGA/Data/New/somatic_binary.csv")


# above tables excluding hypermutators ------------------------------------
count_mutate = somatics_count
counts = c()
mutation = read.csv("hypermutator.csv", skip = 8, nrow = 344)
patients = rownames(somatics_count)
for (i in 1:nrow(mutation)){
  position = grep(mutation$TCGA_Barcode[i],patients)
  counts = append(counts,position)
}
count_mutate = count_mutate[-counts,]
binary_mutate = ifelse(count_mutate>0,1,0)
write.csv(count_mutate, file = "~/TCGA/Data/New/count_mutate.csv")
write.csv(binary_mutate, file = "~/TCGA/Data/New/binary_mutate.csv")
