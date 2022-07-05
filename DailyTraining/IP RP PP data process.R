library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(pheatmap)

test1 <- read.table("GSE151310_CART_IP.UMI_matrix.csv.gz",sep = ",",header = TRUE) #read the dataset
test2 <- read.table("GSE151310_CART_RP.UMI_matrix.csv.gz",sep = ",",header = TRUE)
test3 <- read.table("GSE151310_CART_PP.UMI_matrix.csv.gz",sep = ",",header = TRUE)
genes <- read.table("genes.txt",sep = ",",header = FALSE) %>% unlist


rownames(test1) <- test1[,1]
test1 <- test1[,-1]
rownames(test2) <- test2[,1]
test2 <- test2[,-1]
rownames(test3) <- test3[,1]
test3 <- test3[,-1]

test1[1:5,1:5]

data1 <- CreateSeuratObject(test1,project = "IP")
data2 <- CreateSeuratObject(test2,project = "RP")
data3 <- CreateSeuratObject(test3,project = "PP")


data1[["percent.mt"]] <- PercentageFeatureSet(data1, pattern = "^MT-")
data2[["percent.mt"]] <- PercentageFeatureSet(data2, pattern = "^MT-")
data3[["percent.mt"]] <- PercentageFeatureSet(data3, pattern = "^MT-")#QC with mitochondria genome
VlnPlot(data1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Visualize QC metrics as a violin plot

plot1 <- FeatureScatter(data1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

data1 <- subset(data1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #remove mitochondria leaked died cells
data1 <- NormalizeData(data1, normalization.method = "LogNormalize", scale.factor = 10000)
data2 <- subset(data2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 
data2 <- NormalizeData(data2, normalization.method = "LogNormalize", scale.factor = 10000)
data3 <- subset(data3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 
data3 <- NormalizeData(data3, normalization.method = "LogNormalize", scale.factor = 10000)

re1 <- as.data.frame(data1@assays$RNA@data)
re2 <- as.data.frame(data2@assays$RNA@data)
re3 <- as.data.frame(data3@assays$RNA@data)

ret1 <- as.data.frame(t(re1))
ret2 <- as.data.frame(t(re2))
ret3 <- as.data.frame(t(re3))

plot(ret1$CD8A)
plot(ret3$CD3E)

fin1 <- ret1[which(ret1$CD8A > 0.5 & ret1$CD3E > 0.5),]
fin2 <- ret2[which(ret2$CD8A > 1 & ret2$CD3E > 0.5),]
fin3 <- ret3[which(ret3$CD8A > 1 & ret3$CD3E > 0.5),]

write.csv(fin1[,genes],'D:/RData/IP.csv')
write.csv(fin2[,genes],'D:/RData/RP.csv')
write.csv(fin3[,genes],'D:/RData/PP.csv')


write.csv(ret1,'D:/RData/IP_CART.csv')
write.csv(ret2,'D:/RData/RP_CART.csv')
write.csv(ret3,'D:/RData/PP_CART.csv')


test4 <- read.table("GSE151310_Tcells_PP.UMI_matrix.csv.gz",sep = ",",header = TRUE) #read the dataset
test5 <- read.table("GSE151310_Tcells_RP.UMI_matrix.csv.gz",sep = ",",header = TRUE)

rownames(test4) <- test4[,1]
test4 <- test4[,-1]
rownames(test5) <- test5[,1]
test5 <- test5[,-1]

data4 <- CreateSeuratObject(test4,project = "PP")
data5 <- CreateSeuratObject(test5,project = "RP")

data4[["percent.mt"]] <- PercentageFeatureSet(data4, pattern = "^MT-")
data5[["percent.mt"]] <- PercentageFeatureSet(data5, pattern = "^MT-")

data4 <- subset(data4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #remove mitochondria leaked died cells
data4 <- NormalizeData(data4, normalization.method = "LogNormalize", scale.factor = 10000)
data5 <- subset(data5, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 
data5 <- NormalizeData(data5, normalization.method = "LogNormalize", scale.factor = 10000)


re4 <- as.data.frame(data4@assays$RNA@data)
re5 <- as.data.frame(data5@assays$RNA@data)
ret4 <- as.data.frame(t(re4))
ret5 <- as.data.frame(t(re5))

write.csv(ret4,'D:/RData/PP_Tcell.csv')
write.csv(ret5,'D:/RData/RP_Tcell.csv')