# Title: Multiclass classification for prediction
# Author: Aijin Wang, Contact details: aijinwang3@gmail.com
# Script info:This file contains different approaches to fit lasso logistic 
# regression to predict the cancer type from gene mutations.
#   1. use cv.glmnet on the whole dataset, then predict on the whole dataset
#   2. use cv.glmnet in each fold for training, find lambda, then predict on the test fold
#   3. one vs. rest binary classification
#   4. predict on test fold with all lambdas, pick the best model outside the cv
#  Applications:
#   1. beta coefficients for each gene
#   2. classification analysis
#   3. ROC curve
# 2018 QSURE program at Memorial Sloan Kettering Cancer Center

setwd("~/TCGA/Data/New")
set.seed(1)

library(glmnet)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(ROCR)

type = read.csv("type_h.csv")
type = type[,-1]

x = as.matrix(type[,2:ncol(type)])
y = as.matrix(type[,1])

# 1. Cross validation with cv.glmnet, then predict the whole dataset
   # cv.glmnet to train lasso multinomial logit model ----
### don't run(takes forever)
fit.lasso = cv.glmnet(x,y,
                      family = "multinomial",
                      alpha = 1,
                      nfolds = 3,
                      type.measure = "class")
save(fit.lasso, file = "fit_lasso_logit.rds")

   # use the best lambda in cv to predict values on the whole model ----
fit.lasso <- readRDS("fit_lasso_logit.rds")
best.lambda = fit.lasso$lambda.min
cv.error = fit.lasso$cvm[fit.lasso$lambda == fit.lasso$lambda.min]
model.fit = fit.lasso$glmnet.fit
predictions = predict(model.fit,
                      newx = x, 
                      type = "class",
                      s = best.lambda)

 #* extract coefficients from the model -------------------------------------
best.lambda.index = which(fit.lasso$lambda == fit.lasso$lambda.min)
coefs = model.fit$beta
coefs.df = matrix(NA, 33,299)
  
for (i in 1:length(coefs)){
  cancer.type = coefs[[names(coefs[i])]]
  coefs.df[i,] = as.matrix(cancer.type)[,best.lambda.index]
}

rownames(coefs.df) = names(coefs)
colnames(coefs.df) = rownames(as.matrix(cancer.type))

coefs.df = round(coefs.df,2)


   # order by the number of non-zero coefficients for genes ------------------
binary = ifelse(coefs.df != 0, 1, 0)
ranking = colSums(binary)

dcancer = dist(coefs.df)
hcancer = hclust(dcancer)
coefs.df = coefs.df[hcancer$order, order(ranking)]


   # heatmap of coefficients -------------------------------------------------
meltdata = melt(coefs.df, 
                value.name = "Coefficients", 
                varnames = c("Cancer_Type", "Gene"))
meltdata$Gene = factor(meltdata$Gene, 
                       levels = unique(meltdata$Gene))
meltdata$Cancer_Type = factor(meltdata$Cancer_Type, 
                       levels = unique(meltdata$Cancer_Type))

gene_heat = function(df){
  ggplot(data = df,
         aes(x=Cancer_Type,
             y=Gene,
             fill = Coefficients)) +
    geom_tile() + 
    geom_text(aes(label = Coefficients),size = 2) +
    scale_fill_gradient2(mid = "white",
                         high = "#144182",
                         low= "#8b0000",
                         midpoint = 0,
                         limits = c(-2.776,6.225)) +
    theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8),
          axis.text.y = element_text(size = 5),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none") 
}
ggsave(gene_heat(meltdata),file="~/TCGA/Graphs/heat_coefs.pdf",width=8.5,height = 44)




# 2. Cross validation with cv.glmnet inside 3-fold for loop ----
   # in each fold use cv function to calculate the optimal lambda ------------
n = nrow(type)
k = 3
samp = sample(rep(1:k, length.out = n))
pred.table = matrix(, nrow = 0, ncol = 2)
for (i in 1:k){
  train.x = x[!(samp == k),]
  train.y = as.matrix(y[!(samp == k),])
  test.x = x[samp == k,]
  test.y = as.matrix(y[samp == k,])
  fit.cv = cv.glmnet(train.x, train.y,
                     family = "multinomial",
                     alpha = 1,
                     nfolds = 3,
                     type.measure = "class")
  fit = fit.cv$glmnet.fit
  best.lambda = fit.cv$lambda.min
  pred = predict(fit,
                 newx = test.x,
                 type = "class",
                 s = best.lambda)
  pred.table = rbind(pred.table, data.frame(test.y, pred))
}
pred.table[] = lapply(pred.table, as.character)
write.csv(pred.table, file = "cv_cv.csv", row.names = FALSE)



 # *classification analysis----

   # sensitivity -------------------------------------------------------------
classification_test = function(cancer){
  TP = sum(pred.table$test.y == cancer & pred.table$X1 == cancer)
  FP = sum(pred.table$test.y != cancer & pred.table$X1 == cancer)
  FN = sum(pred.table$test.y == cancer & pred.table$X1 != cancer)
  TN = sum(pred.table$test.y != cancer & pred.table$X1 != cancer)
  pos_pred = TP / (TP + FP)
  neg_pred = TN / (FN + TN)
  sensitivity = TP / (TP + FN)
  specificity = TN / (FP + TN)
  ACC = (TP + TN) / (TP + FN + TN + FP)
  ER = (FN + FP) / (TP + FN + TN + FP)
  return(c(pos_pred, neg_pred, 
           sensitivity, specificity,
           ACC, ER))
}

cancer_names = unique(pred.table$test.y)
classification.df = data.frame()
for (i in 1:length(cancer_names)){
  classification.df = rbind(classification.df, 
                            classification_test(cancer_names[i]))
}
colnames(classification.df) = c("Positive Predictive Value", 
                                "Negative Predictive Value", 
                                "Sensitivity", "Specificity",
                                "Accuracy", "Error Rate")
classification.df$cancer = cancer_names


   # plot classification result ----------------------------------------------
plotting = function(column){
  ggplot(classification.df, 
         aes(x = reorder(cancer, -get(column)), 
             y = get(column))) + 
    geom_point() +
    labs(x = "Cancer Type",
         y = get("column")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

grid.arrange(plotting("Positive Predictive Value"), 
             plotting("Negative Predictive Value"), 
             plotting("Sensitivity"), 
             plotting("Specificity"),
             plotting("Accuracy"),
             plotting("Error Rate"))


# 3. One vs. rest binary classification ----
cancer_names = unique(type$X)
n = nrow(type)
k = 3
samp = sample(rep(1:k, length.out = n))

result = list()
for (type.id in 1:33){
  cancer.predict.prob = c()
  cancer.predict.class = c()
  cancer.actual = c()
  for (k in 1:3){
    traind = type[!(samp == k),]
    testd = type[samp == k,]

    cancer = as.factor(traind$X == cancer_names[type.id])
    lasso.model = try(cv.glmnet(as.matrix(traind[,-1]), as.matrix(cancer),
                          family = "binomial",
                          nfolds = 3,
                          type.measure = "class"))
    lasso.predict = try(predict(lasso.model,
                          newx = as.matrix(testd[,-1]),
                          s = lasso.model$lambda.min,
                          type = "response"))
    lasso.predict.class = try(predict(lasso.model,
                                newx = as.matrix(testd[,-1]),
                                s = lasso.model$lambda.min,
                                type = "class"))
    actual.class = try(testd$X == cancer_names[type.id])
    cancer.predict.prob = try(append(cancer.predict.prob, lasso.predict))
    cancer.predict.class = try(append(cancer.predict.class, lasso.predict.class))
    cancer.actual = try(append(cancer.actual, actual.class))
  }
  result[[type.id]] = try(list(cancer.actual, 
                                      cancer.predict.class,
                                      cancer.predict.prob))
  result[[type.id]][[2]] = ifelse(result[[type.id]][[2]] == "TRUE",TRUE, FALSE)
}
# saveRDS(result, file = "cvresult_v2.rds")

 #* ROC ----
aucs = c()
plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),
     ylab='True Positive Rate',
     xlab='False Positive Rate',
     bty='n')

for(type.id in 1:33){
  cancer.result = result[[type.id]]
  pred = prediction(cancer.result[[3]], cancer.result[[1]])
  perf = performance(pred, "tpr", "fpr")
  roc.x = unlist(perf@x.values)
  roc.y = unlist(perf@y.values)

  lines(roc.y ~ roc.x, col = type.id + 1, lwd = 2)

  nbauc = performance(pred, "auc")
  nbauc = unlist(slot(nbauc, "y.values"))
  aucs[type.id] = nbauc
}

# plot AUC
plot.data = data.frame(cancer_names, aucs)

ggplot(plot.data, 
       aes(x = reorder(cancer_names, -aucs), 
           y = aucs)) + 
  geom_point() +
  labs(x = "Cancer Type",
       y = "AUC") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

   # classification results and plots ----
classification_test_v2 = function(cancer){
  TP = sum(result[[cancer]][[1]] == TRUE & result[[cancer]][[2]] == TRUE)
  FP = sum(result[[cancer]][[1]] != TRUE & result[[cancer]][[2]] == TRUE)
  FN = sum(result[[cancer]][[1]] == TRUE & result[[cancer]][[2]] != TRUE)
  TN = sum(result[[cancer]][[1]] != TRUE & result[[cancer]][[2]] != TRUE)
  pos_pred = TP / (TP + FP)
  neg_pred = TN / (FN + TN)
  sensitivity = TP / (TP + FN)
  specificity = TN / (FP + TN)
  ACC = (TP + TN) / (TP + FN + TN + FP)
  ER = (FN + FP) / (TP + FN + TN + FP)
  ranking = pos_pred + sensitivity - 1
  return(c(pos_pred, neg_pred, 
           sensitivity, specificity,
           ACC, ER, ranking))
}

classification.df.v2 = data.frame()
for (i in 1:length(cancer_names)){
  classification.df.v2 = rbind(classification.df.v2, 
                            classification_test_v2(i))
}
colnames(classification.df.v2) = c("Positive Predictive Value", 
                                "Negative Predictive Value", 
                                "Sensitivity", "Specificity",
                                "Accuracy", "Error Rate", "Ranking")
classification.df.v2$cancer = cancer_names

plotting_v2 = function(column){
  ggplot(classification.df.v2, 
         aes(x = reorder(cancer, -get(column)), 
             y = get(column))) + 
    geom_point() +
    labs(x = "Cancer Type",
         y = get("column")) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

grid.arrange(plotting_v2("Positive Predictive Value"), 
             plotting_v2("Negative Predictive Value"), 
             plotting_v2("Sensitivity"), 
             plotting_v2("Specificity"),
             plotting_v2("Accuracy"),
             plotting_v2("Error Rate"))

plotting_v2("Ranking")

# 4. Cross validation with different assigned lambdas ----
k = 3
samp = sample(rep(1:k, length.out = n))
pred.table = matrix(, nrow = 0, ncol = 50)
real.y = c()
for (i in 1:k){
  
  train.x = x[!(samp == k),]
  train.y = as.matrix(y[!(samp == k),])
  
  test.x = x[samp == k,]
  test.y = as.matrix(y[samp == k,])
  
  fit = glmnet(train.x, train.y,
               family = "multinomial",
               alpha = 1,
               lambda = seq(1.92*10^(-5), 0.1753, length = 50))
  pred = predict(fit,
                 newx = test.x,
                 type = "class")
  pred = as.data.frame(pred)
  pred.table = bind_rows(pred.table, pred)
  real.y = append(real.y, test.y)
}
pred.table = pred.table[-1,]
#write.csv(pred.table, file = "cv.csv", row.names = FALSE)

colnames(pred.table) = sort(seq(1.92*10^(-5), 0.1753, length = 50),
                            decreasing = TRUE)[1:44]
pred.table$real.y = real.y

