# Title: Frequency Heatmap
# Author: Aijin Wang, Contact details: aijinwang3@gmail.com
# Script info: Heatmap of mutation freq rate for each gene in each cancer type 
#              1. alphabetical order
#              2. hierarchiel clustering order
#              3. pathway gene alphabetical order
# 2018 QSURE program at Memorial Sloan Kettering Cancer Center

setwd("~/TCGA/Data/Raw")
library(readxl)
library(plyr)
library(tidyverse)
library(gplots)
library(reshape2)

# map the patient ID with tumor type in binary dataset -----
mapping = read_excel("sampleID.xlsx")
somatics_binary = read.csv("~/TCGA/Data/New/binary_mutate.csv")
patients = somatics_binary[,1]
somatics_binary[,1] = as.character(somatics_binary[,1])
for (i in 1:nrow(mapping)){
  position = grep(mapping$ParticipantBarcode[i], patients)
  if (!identical(position, integer(0))){
    somatics_binary[position,1] = mapping$Study[i]
  }
}


# delete rows with unknown tumor types -----
somatics_binary = somatics_binary[nchar(somatics_binary[,1]) < 10,]
write.csv(somatics_binary, file = "~/TCGA/Data/New/type_h.csv")


# calculate ratio of mutated genes in each tumor type -----
ratio = function(gene){
  new = data.frame(tumor = somatics_binary[,1],
                   gene = somatics_binary[,gene])
  gene_table = table(new)
  gene_table = cbind(gene_table, 
                     ratio = gene_table[,"1"] / rowSums(gene_table)*100)
  return(gene_table[,"ratio"])
}

ratio_table = c()
for (i in 2:ncol(somatics_binary)){
  ratio_table = cbind(ratio_table, ratio(colnames(somatics_binary)[i]))
}
colnames(ratio_table) = colnames(somatics_binary)[2:ncol(somatics_binary)]
write.csv(ratio_table, file = "~/TCGA/Data/New/ratio_table.csv")


# mutation frequency ratio heatmap -----
ratio_table = round(ratio_table, 1)

gene_heat = function(df){
  ggplot(data = df,
         aes(x=Cancer_Type,
             y=Gene,
             fill = Proportion)) +
    geom_tile() + 
    geom_text(aes(label = Proportion),
              size = 2) +
    scale_fill_gradient2(low = "white",
                         mid = "steelblue",
                         high="#2f5c33",
                         midpoint = 50) +
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8),
          axis.text.y = element_text(size = 5),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") 
}


# 1.alphabetical order
meltdata = melt(ratio_table, 
                value.name = "Proportion", 
                varnames = c("Cancer_Type", "Gene"))

ggsave(gene_heat(meltdata),
       file="~/TCGA/Graphs/heat_original_mutate.pdf",
       width=8.5,
       height = 44)


# 2.hierarchiel clustering order
dcancer = dist(ratio_table)
hcancer = hclust(dcancer)
dgene = dist(t(ratio_table))
hgene = hclust(dgene)

cluster_ratio = ratio_table
cluster_ratio = cluster_ratio[hcancer$order, hgene$order]
melt_cluster = melt(cluster_ratio, 
                    value.name = "Proportion",
                    varnames = c("Cancer_Type", "Gene"))
melt_cluster$Cancer_Type = factor(melt_cluster$Cancer_Type, 
                                  levels = unique(melt_cluster$Cancer_Type))
melt_cluster$Gene = factor(melt_cluster$Gene, 
                           levels = unique(melt_cluster$Gene))

ggsave(gene_heat(melt_cluster),
       file="~/TCGA/Graphs/heat_cluster_mutate.pdf",
       width=8.5,
       height = 44)


# 3.pathway gene
pathways = read.csv("pathway.csv")
pathway = unique(pathways$Gene)
pathway_ratio = ratio_table[,colnames(ratio_table) %in% pathway]

dcancer_p = dist(pathway_ratio)
hcancer_p = hclust(dcancer_p)

pathway_ratio = pathway_ratio[hcancer_p$order,]
melt_p = melt(pathway_ratio, 
              value.name = "Proportion",
              varnames = c("Cancer_Type", "Gene"))
melt_p$Cancer_Type = factor(melt_p$Cancer_Type, 
                            levels = unique(melt_p$Cancer_Type))
melt_p$Gene = factor(melt_p$Gene, 
                     levels = unique(melt_p$Gene))

ggsave(gene_heat(melt_p),
       file="~/TCGA/Graphs/heat_cluster_pathway.pdf",
       width=8.5,
       height = 44)
