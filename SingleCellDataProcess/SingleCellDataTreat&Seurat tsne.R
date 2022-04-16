library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(pheatmap)

patient <- read.table("GSE120575_patient_ID_single_cells.txt",header = F)
header <- read.table("GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt",sep = "\t",header = FALSE,nrows = 2)
test1 <- read.table("GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt",sep = "\t",header = FALSE,skip = 2) #read the dataset
test2 <- read.table("GSE151310_CART_RP.UMI_matrix.csv.gz",sep = ",",header = TRUE)
test3 <- read.table("GSE151310_CART_IP.UMI_matrix.csv.gz",sep = ",",header = TRUE)
genes <- read.table("genes.txt",sep = ",",header = FALSE) %>% unlist

rownames(test1) <- test1[,1]
test1 <- test1[,-1]
rownames(test2) <- test2[,1]
test2 <- test2[,-1]
rownames(test3) <- test3[,1]
test3 <- test3[,-1]

data1 <- test1[,which(unlist(test1["MYLIP",])!=0)]
data2 <- test2[,which(unlist(test2["MYLIP",])!=0)]
data3 <- test3[,which(unlist(test3["MYLIP",])!=0)]

dataheatmap <- cbind(data1,data2,data3)
pheatmap(data1[c("IFNG","MYLIP","GZMB"),],scale = "row")


ggplot(data4,aes(`IFNG`,`MYLIP`))+geom_dotplot()

data1 <- CreateSeuratObject(counts = data1,project = "PP",min.cells = 3, min.features = 200) #select gene express in more than 3 cells and cells express more than 200 genes
data2 <- CreateSeuratObject(counts = data2,project = "RP",min.cells = 3, min.features = 200)
data3 <- CreateSeuratObject(counts = data3,project = "IP",min.cells = 3, min.features = 200)
data <- merge(data1, y = data2, add.cell.ids = c("PP", "RP"), project = "CART")

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


saveRDS(data, file = "../pbmc_tutorial.rds")

pdf("",height = 4,width = 6)
plot1 + plot2
dev.off()
