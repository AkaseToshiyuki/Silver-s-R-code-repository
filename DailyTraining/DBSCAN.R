library(RTCGA)
library(RTCGA.rnaseq)
library(ggplot2)
library(stringr)
library(factoextra)
library(dbscan)
library(fpc)

exp <- expressionsTCGA(BRCA.rnaseq) %>% as.data.frame()
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp <- exp[,-20532]
exp <- scale(exp)

expr <- exp[,which(colSums(exp)!=0)]


kNNdistplot(expr,k = 10)

DB <- dbscan(expr,eps = 200,MinPts = 10,scale = F)
print(DB)
fviz_cluster(DB,expr,repel = T)
