library(Seurat)
library(ggplot2)
library(Rtsne)
library(dplyr)
#载入包
data <- read.csv("Data.csv")
#读取数据框
dim(data)
#确定数据格式，为基因数（维度）：样品数
#如果不是就转置为以上格式
data1 <- t(data)
data <- data.frame(data1)
dim(data)
#转置
cip1 <-CreateSeuratObject(data)
cip1
#应该如下
#An object of class Seurat 
#6 features across 4751 samples within 1 assay 
#Active assay: RNA (6 features, 0 variable features)
all.genes <- rownames(x = cip1)
#给类命名
cip2<- ScaleData(object =cip1, features = all.genes)
#数据缩放
FindVariableFeatures(cip2)
cip2<- RunPCA(object = cip2, features = all.genes)
#PCA分析
VizDimLoadings(object = cip2, dims = 1:2, reduction = "pca")
#梯度图生成
DimHeatmap(object = cip2, dims = 1:5, cells = 100, balanced = TRUE)
#热图生成
cip2 <- JackStraw(object = cip2, num.replicate = 100)
cip2 <- ScoreJackStraw(object = cip2, dims = 1:5)
JackStrawPlot(object = cip2, dims = 1:5)
#选择维度，所有的dims = 1:X 中的X都用基因数-1来表示，大于20则统一为20
cip2 <- FindNeighbors(object = cip2, dims = 1:5)
cip2 <- FindClusters(object = cip2, resolution = 0.5)
head(x = Idents(object = cip2), 5)
#聚类分析
cip2 <- RunTSNE(object = cip2, dims = 1:5)
TSNEPlot(object = cip2)
FeaturePlot(cip2, all.genes, pt.size = 1,ncol = 2)
#作图
#以下为手动注释步骤，如果有参考数据集可以不用手动注释，自行寻找解决方案
Monocyte = c("CD4","CD45","CD11b")
DC = c("CD4","CD11b")
Treg = c("CD4","CD25","CD45")
Tfh = c("CD4")
CTL = c("CD8")
Neutrophil = c("CD45","CD11b","Ly6G")
NK = c("CD11b")
subtypes = c(Monocyte,DC,Treg,Tfh,CTL,Neutrophil,NK)

