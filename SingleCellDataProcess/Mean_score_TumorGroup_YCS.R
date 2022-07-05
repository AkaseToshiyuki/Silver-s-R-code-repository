library(dplyr)
library(stringr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(clusterProfiler) 
library(ggplot2)

biosyngene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
uptakegene <- c("LDLR", "NPC1", "NPC2")
SREBF2genes <- c("SREBF2","HMGCR","SQLE")


raw <- read.table("GSE140228_read_counts_Smartseq2.csv.gz",sep = ",",header = T)
colnames(raw) <- str_replace_all(colnames(raw),"\\.","-")
cell <- read.table("GSE140228_cell_info_Smartseq2.tsv.gz",sep = "\t",header = T)
rownames(raw) <- raw[,1]
raw <- raw[,-1]
table(cell$Tissue)
table(cell$celltype_sub)

tumorID <- cell[which(cell$Tissue == "Tumor"),]
table(tumorID$Sample)
table(tumorID$celltype_sub)

tumorCD4 <- tumorID[which(str_sub(tumorID$celltype_sub,1,3)== "CD4"),]
tumorCD8 <- tumorID[which(str_sub(tumorID$celltype_sub,1,3)== "CD8"),]
tumorDC <- tumorID[which(str_sub(tumorID$celltype_sub,1,2)== "DC"),]
tumorB <- tumorID[which(str_sub(tumorID$celltype_sub,1,3)== "Lym"),]
tumorMast <- tumorID[which(str_sub(tumorID$celltype_sub,1,4)== "Mast"),]
tumorMono <- tumorID[which(str_sub(tumorID$celltype_sub,1,4)== "Mono"),]
tumorM <- tumorID[which(str_sub(tumorID$celltype_sub,1,2)== "M蠁"),]
tumorNK <- tumorID[which(str_sub(tumorID$celltype_sub,1,2)== "NK"),]

dttumorCD4 <- raw[,which(colnames(raw) %in% tumorCD4$Barcode)]
dttumorCD8 <- raw[,which(colnames(raw) %in% tumorCD8$Barcode)]
dttumorDC <- raw[,which(colnames(raw) %in% tumorDC$Barcode)]
dttumorB <- raw[,which(colnames(raw) %in% tumorB$Barcode)]
dttumorMast <- raw[,which(colnames(raw) %in% tumorMast$Barcode)]
dttumorMono <- raw[,which(colnames(raw) %in% tumorMono$Barcode)]
dttumorM <- raw[,which(colnames(raw) %in% tumorM$Barcode)]
dttumorNK <- raw[,which(colnames(raw) %in% tumorNK$Barcode)]

tcd4 <- CreateSeuratObject(dttumorCD4,project = "CD4",min.cells = 1, min.features = 20)
tcd8 <- CreateSeuratObject(dttumorCD8,project = "CD8",min.cells = 1, min.features = 20)
tDC <- CreateSeuratObject(dttumorDC,project = "DC",min.cells = 1, min.features = 20)
tB <- CreateSeuratObject(dttumorB,project = "B",min.cells = 1, min.features = 20)
tMast <- CreateSeuratObject(dttumorMast,project = "Mast",min.cells = 1, min.features = 20)
tMono <- CreateSeuratObject(dttumorMono,project = "Mono",min.cells = 1, min.features = 20)
tM <- CreateSeuratObject(dttumorM,project = "M",min.cells = 1, min.features = 20)
tNK <- CreateSeuratObject(dttumorNK,project = "NK",min.cells = 1, min.features = 20)

cd8data <- merge(x = tcd8,y =c(tcd4,tDC,tB,tMast,tMono,tM,tNK),add.cell.ids = c("CD8","CD4","DC","B","Mast","Mono","M","NK"))

cd8data[["percent.mt"]] <- PercentageFeatureSet(cd8data, pattern = "^MT-") 
VlnPlot(cd8data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(cd8data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cd8data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
cd8data <- NormalizeData(cd8data, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(cd8data)
cd8data <- ScaleData(cd8data, features = biosyngene)

re <- as.data.frame(cd8data@assays$RNA@data)
re[1:5,1:5]

metadata <- re[syngene,]
metadata[nrow(metadata)+1,] <- colMeans(metadata)
plotdata <- metadata[nrow(metadata),]
plotdata <- t(plotdata) %>% as.data.frame()
plotdata$Group <- str_extract(rownames(plotdata),".*(?=_[:upper:])")
colnames(plotdata)[1] <- "Score"

ggplot(plotdata,aes(`Group`,`Score`))+ 
  geom_violin(aes(fill = `Group`))+
  geom_boxplot(width = 0.4,fill = "white",outlier.alpha = 0)+
  geom_jitter(size = 0.2, width = 0.1,shape = 16,alpha = 0.5,color = "grey20")+
  theme_classic()
