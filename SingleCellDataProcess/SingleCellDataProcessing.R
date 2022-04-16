library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(stringr)
library(ggsignif)


patient <- read.table("GSE120575_patient_ID_single_cells.txt",header = F)
header <- read.table("GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt",sep = "\t",header = FALSE,nrows = 2)
test1 <- read.table("GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt",sep = "\t",header = FALSE,skip = 2) #read the dataset
genes <- read.table("genes.txt",sep = ",",header = FALSE) %>% unlist

colnames(test1) <- header[1,] 

responder <- patient[which(patient$V10=="Responder"&str_detect(patient$V9,"Post")),]
grouplist <- as.data.frame(aggregate(responder$V9,list(responder$V9),length)) #or use table($) to get the list of post_
grouplist <- unlist(grouplist[,-2])

colnames(header) <- header[2,]
header <- header[-2,]
candidate <- header[,colnames(header) %in% grouplist] 
cellid<- unlist(candidate[1,])


header <- read.table("GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt",sep = "\t",header = FALSE,nrows = 2)
nonresponder <- patient[which(patient$V10=="Non-responder"&str_detect(patient$V9,"Post")),]
nongrouplist <- as.data.frame(aggregate(nonresponder$V9,list(nonresponder$V9),length))
nongrouplist <- unlist(nongrouplist[,-2])

noncandidate <- header[,colnames(header) %in% nongrouplist] 
noncellid<- unlist(noncandidate[1,])



datar <- test1[,colnames(test1) %in% cellid]
rownames(datar) <- test1[,1]
datar[1:5,1:5]
pheatmap(datar[c("IFNG","MYLIP","GZMB"),],scale = "row")
datag <- as.data.frame(t(datar))
datag[1:5,1:5]
plot(datag$CD8A)
datar8 <- datag[which(datag$CD8A > 8),]
datar8[1:5,1:5]
plot(datar8$CD8A)



datanr <- test1[,colnames(test1) %in% noncellid]
rownames(datanr) <- test1[,1]
datanr[1:5,1:5]
pheatmap(datanr[c("IFNG","MYLIP","GZMB"),],scale = "row")
datang <- as.data.frame(t(datanr))
datang[1:5,1:5]
plot(datang$CD8A)
datanr8 <- datang[which(datang$CD8A > 8),]
datanr8[1:5,1:5]
plot(datar8$CD8A)


rcd8 <- as.data.frame(datar8[,"MYLIP"])
rcd8 <- rcd8[rcd8$`datar8[, "MYLIP"]`!=0,]
nrcd8 <- as.data.frame(datanr8[,"MYLIP"])
nrcd8 <- nrcd8[nrcd8$`datanr8[, "MYLIP"]`!=0,]

write.csv(rcd8,'D:/RData/responder_CD8+_MYLIP.csv')
write.csv(nrcd8,'D:/RData/nonresponder_CD8+_MYLIP.csv')


threedatar <- datar8[,genes]
threedatanr <- datanr8[,genes]

write.csv(threedatar,'D:/RData/responder_CD8+.csv')
write.csv(threedatanr,'D:/RData/nonresponder_CD8+.csv')

ggplot(datar8,mapping = aes(`IFNG`,`MYLIP`)) + geom_point()







#dont use follow code

namereplace <- colnames(datar)
namereplace <- str_replace_all(namereplace,"_","")
colnames(datar) <- namereplace
datar[1:5,1:5]



genereplace <- rownames(datar)
genereplace <- str_replace_all(genereplace,"_","")
rownames(datar) <- genereplace
datar[1:5,1:5]





data <-CreateSeuratObject(counts = datar,project = "ResponderPost",min.cells = 5, min.features = 200)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-") #QC with mitochondria genome
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Visualize QC metrics as a violin plot

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #remove mitochondria leaked died cells
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000) #multi 10000 and get log to normalize

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000) #identifiy of features
top30 <- head(VariableFeatures(data), 30)
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = genes, repel = TRUE)
plot2

all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data, features = VariableFeatures(object = data))
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca")
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)

data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:20)
JackStrawPlot(data, dims = 1:15)

data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.5)
data <- RunTSNE(data, dims = 1:10)

DimPlot(data, reduction = "tsne")
DimPlot(data, reduction = "tsne", group.by = "orig.ident")
FeaturePlot(data, features = c( "SREBF2", "HMGCR", "SQLE", "LDLR"))
FeaturePlot(data, features = c( "NR1H3", "NR1H2", "ABCA1", "ABCG1","MYLIP"))
FeaturePlot(data, features = c("IFNG","ABCG1","GZMB"))
