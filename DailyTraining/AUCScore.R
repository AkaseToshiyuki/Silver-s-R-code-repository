library(dplyr)
library(stringr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(AUCell) 
library(clusterProfiler) 
library(ggplot2)


genes <- read.table("Cholesterol metabolism genes_2019-05-30.csv",sep = ",",header = T,col.names = c("Gene","Pathway"))
biosyngene <- genes[which(genes$Pathway == "Biosynthesis"),1]
uptakegene <- c("LDLR","SCARB1")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2")


raw <- read.table("GSE140228_read_counts_Smartseq2.csv.gz",sep = ",",header = T)
colnames(raw) <- str_replace_all(colnames(raw),"\\.","-")
cell <- read.table("GSE140228_cell_info_Smartseq2.tsv.gz",sep = "\t",header = T)
rownames(raw) <- raw[,1]
raw <- raw[,-1]
table(cell$Tissue)

tumorID <- cell[which(cell$Tissue == "Tumor"),]
normalID <- cell[which(cell$Tissue == "Normal"),]
table(tumorID$Sample)
table(normalID$Sample)
table(tumorID$celltype_sub)
table(normalID$celltype_sub)

tumorCD8 <- tumorID[which(str_sub(tumorID$celltype_sub,1,3)== "CD8"),]
normalCD8 <- normalID[which(str_sub(normalID$celltype_sub,1,3)== "CD8"),]
dttumorCD8 <- raw[,which(colnames(raw) %in% tumorCD8$Barcode)]
dtnormalCD8 <- raw[,which(colnames(raw) %in% normalCD8$Barcode)]

tcd8 <- CreateSeuratObject(dttumorCD8,project = "tumor",min.cells = 1, min.features = 20)
ncd8 <- CreateSeuratObject(dtnormalCD8,project = "normal",min.cells = 1, min.features = 20)
cd8data <- merge(x = tcd8,y = ncd8,add.cell.ids = c("tumor","normal"))
# Mulitigroup
# cd8data <- merge(groupA,y = c(groupB,groupC,groupD),add.cell.ids = c("groupA","groupB","groupC","groupD"))


cd8data[["percent.mt"]] <- PercentageFeatureSet(cd8data, pattern = "^MT-") 
VlnPlot(cd8data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(cd8data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cd8data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
cd8data <- NormalizeData(cd8data, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(cd8data)
cd8data <- ScaleData(cd8data, features = all.genes)
#fpkm



re <- as.matrix(cd8data@assays$RNA@scale.data)
re[1:5,1:5]

cells_rankings <- AUCell_buildRankings(re) ##基因排序
cells_AUC <- AUCell_calcAUC(LXRgene, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05) ##计算AUC值。
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1,assign=TRUE) ##挑选阈值

dt <- as.data.frame(cells_AUC@assays@data@listData)
plotdata <- t(dt) %>% as.data.frame()
# <- rbind(t(normalAUC),t(tumorAUC)) %>% as.data.frame()
plotdata$Group <- str_extract(rownames(plotdata),"(?<=\\.).*(?=_[:upper:])")


ggplot(plotdata,aes(`Group`,`geneSet`)) + geom_violin(aes(fill = `Group`))+geom_boxplot(width = 0.05,fill = "white",outlier.alpha = 0)+
  theme_classic()
  
