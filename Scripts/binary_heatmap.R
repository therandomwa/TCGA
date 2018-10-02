# Title: Binary Mutation Heatmap
# Author: Aijin Wang, Contact details: aijinwang3@gmail.com
# Script info: binary heatmap of hierchical clustering at individual gene level
# 2018 QSURE program at Memorial Sloan Kettering Cancer Center

setwd("~/TCGA/Data/New")
set.seed(1)
library(pheatmap)


type = read.csv("type_h.csv")
type = type[,-1]
type = as.matrix(type)
type = t(type)


# exclude hypermutators -----
counts = type[2:nrow(type),]
class(counts) = "numeric"
sums = colSums(counts)
clusters = kmeans(sums,2)$cluster
table(clusters)
type = type[,clusters == 1]
cancer = type[1,]
colnames(type) = paste("Patient",1:ncol(type), sep = "")

# create annotation table -----
categories = data.frame(cancer_type = cancer)
rownames(categories) = colnames(type)
type = type[-1,]
class(type) = "numeric"
dim(type)


# specify colors for legends -----
set.seed(11)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
my_color = sample(color,33)
names(my_color) = levels(categories$cancer_type)
cancer_col = list(cancer_type = my_color)


# heatmap -----
pheatmap(type, 
         color = c("white","black"),
         annotation = categories,
         annotation_colors = cancer_col,
         cellwidth = 0.15,
         cellheight = 8,
         treeheight_row = 0,
         treeheight_col = 0,
         legend = F,
         clustering_distance_rows = "binary",
         clustering_distance_cols = "binary",
         clustering_method = "ward.D2",
         show_colnames = F,
         fontsize = 6.5)
