# Title: Scatterpie/beachball
# Author: Aijin Wang, Contact details: aijinwang3@gmail.com
# Script info: Beachball(scatterpie) graph to analyze the relationship between 
# within cluster homogeity and dominant cancer proportions.
# 2018 QSURE program at Memorial Sloan Kettering Cancer Center

setwd("~/TCGA/Data/New")

library(tidyverse)
library(ggforce)
library(scatterpie)
library(cluster)
library(dplyr)

set.seed(1)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
my_color = sample(color,33)

type = read.csv("type_h.csv")
type = type[,-1]
type = as.matrix(type)
rownames(type) = type[,1]
type = type[,-1]

# gather silouette and cluster information ----
d = dist(type, method = "binary")
h = hclust(d, method = "ward.D2")
groups = cutree(h, k = 33)

s = silhouette(groups, d)
clusters = data.frame(Cancer = rownames(type),
                      Groups = s[,"cluster"],
                      Sil_value = s[,"sil_width"])


# calculate the cancer proportion for each pie chart ----
pie_prop = function(x){
  counts = as.data.frame(table(x$Cancer))
  return (counts$Freq)
}
pie_data = data.frame(daply(clusters, .(Groups), pie_prop))
colnames(pie_data) = levels(clusters$Cancer)
pie_data$sums = rowSums(pie_data)/sum(pie_data)
pie_data$sums = pie_data$sums^(1/3)/9


# calculate the mean silouette value for each cluster ----
mean_sil = daply(clusters, .(Groups),function(x){mean(x$Sil_value)})


# find the most dominant cancer type in each cluster ----
dominant_prop = function(x){
  counts = as.data.frame(table(x$Cancer))
  counts$Prop = counts$Freq / sum(counts$Freq)
  props = counts$Prop[which(counts$Prop == max(counts$Prop))]
  name = counts$Var1[which(counts$Prop == max(counts$Prop))]
  props = props[1]
  name = as.character(name[1])
  return(c(name,props))
}
max_data = daply(clusters, .(Groups), dominant_prop)
max_prop = as.numeric(max_data[,2])
max_name = max_data[,1]


# combine all the data ----
axis_data = data.frame(mean_sil = mean_sil, 
                       max_prop = max_prop,
                       max_name = max_name)
axis_data$clust = factor(1:33)
axis_data = cbind(axis_data, pie_data)


# convert the data type ----
for (i in 4:ncol(axis_data)){
  axis_data[,i] = as.numeric(axis_data[,i])
}


# scale using each cluster's proportion ----
ggplot() + 
  geom_scatterpie(data = axis_data,
                  aes(x = mean_sil,
                      y = max_prop, 
                      group = clust,
                      r = sums),
                  cols = levels(clusters$Cancer),
                  size = 0.3) +
  geom_text(data = axis_data,
            aes(x = mean_sil - 0.05,
                y = max_prop + 0.04,
                label = max_name),
            size = 2.5) +
  coord_equal() +
  scale_fill_manual(values = my_color) +
  labs(x = "Mean Silhouette Width",
       y = "Dominant Cancer Type Proportion") +
  theme_bw() +
  theme(legend.text = element_text(size = 7),
        legend.key.size = unit(0.8,"line")) +
  guides(fill = guide_legend(ncol = 1))


# pick mixed group to analyze ----

# mixed clusters number
small_prop = axis_data[(axis_data$max_prop < 0.25) 
                       & axis_data$mean_sil < 1,]$clust
# samples that are in mixed clusters
small_data = type[clusters$Groups %in% small_prop,]  %>% as.data.frame()
Group = as.numeric(clusters[clusters$Groups %in% small_prop, "Groups"])
small_data = cbind(small_data, Group)
for (i in 1:ncol(small_data)){
  small_data[,i] = as.numeric(as.character(small_data[,i]))
}
# count mutation rate for each gene
small_count = daply(small_data, .(Group), function(x){colSums(x)})
small_count = small_count[,-300]



# pick a mixed cluster  ----
n = 1
m = as.numeric(rownames(small_count)[n])
# count histogram
name = data.frame(names = colnames(small_count), value = small_count[n,])
name = name[order(name$value),]
name$names = factor(name$names, levels = unique(name$names))
ggplot(data = name) + 
  geom_bar(aes(x = names, y = value), stat = "identity") +
  labs(x = "Gene",
       y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))

# zoomed in pie chart
pien_data = axis_data[m, 5:(ncol(axis_data)-2)]
pien_df = data.frame(cancer = colnames(pien_data), 
                     count = as.numeric(pien_data[1,]))
maxn_name = axis_data[m,]$max_name
ggplot(data = pien_df, aes(x = "", y = count, fill = cancer)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = my_color) +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "none")
