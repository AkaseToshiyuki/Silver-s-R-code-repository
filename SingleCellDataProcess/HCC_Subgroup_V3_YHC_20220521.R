library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)

SREBPgene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
               "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
               "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24","LDLR")
acquisitiongene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                     "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                     "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24","LDLR", "NPC1", "NPC2")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
uptakegene <- c("SCARB1","LDLR", "NPC1", "NPC2")
biosyngene <- c("ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24")
CDmarker <- c("CD8A","CD8B","CD4","CD3E")
SREBF2genes <- c("SREBF2","HMGCR","SQLE","LDLR")
functiongene <- c("IFNG", "IL10", "GZMA", "GZMB", "PRF1","TNFSF10", "FASLG", "TNF", "IL2","CD44", "IL2RA")
proliferation <- read.table("ProliferationHCC_gene.csv",header = T,sep = ",")
proliferation <- proliferation[,2]

cell <- read.table("GSE140228_cell_info_Smartseq2.tsv.gz",sep = "\t",header = T)
load("TPMData_HCC.RData")
raw <- t(TPM) %>% as.data.frame()
raw[1:5,1:5]
rm(TPM)
gc()

cell <- cell[which(cell$Tissue == "Tumor"),]
table(cell$celltype_sub)

CD4 <- cell[which(str_detect(cell$celltype_sub,"CD4-")),]
CD8 <- cell[which(str_detect(cell$celltype_sub,"CD8-")),]
DC <- cell[which(str_detect(cell$celltype_sub,"DC-")),]
B <- cell[which(str_detect(cell$celltype_sub,"Lymphoid-B")),]
M <- cell[which(str_detect(cell$celltype_sub,"M蠁-")),]
Mono <- cell[which(str_detect(cell$celltype_sub,"Mono-")),]
NK <- cell[which(str_detect(cell$celltype_sub,"NK-")),]
Mast <- cell[which(str_detect(cell$celltype_sub,"Mast-")),]

CD4$Group <- "CD4 T Cell"
CD8$Group <- "CD8 T Cell"
DC$Group <- "DC"
B$Group <- "B cell"
M$Group <- "Macrophages"
Mono$Group <- "Monocytes"
NK$Group <- "NK Cell"
Mast$Group <- "Mast Cell"

cellG <- bind_rows(CD4,CD8,DC,B,M,Mono,NK,Mast)

raw <- raw[which(rownames(raw) %in% cellG$Barcode),]
rownames(cellG) <- cellG$Barcode
cellG <- cellG[rownames(raw),]

geneused <- c(uptakegene,biosyngene)
raw <- raw[,geneused]
scaledata <- as.data.frame(scale(raw))
scaledata$Group <- cellG$Group

scaledata$Uscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% uptakegene)])
scaledata$Sscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% biosyngene)])

plotdata <- scaledata[,c("Uscore","Sscore","Group")]

ggplot(plotdata,aes(`Group`,`Sscore`)) +
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic() +
  stat_compare_means(method = "anova",label = "p.format")

ggplot(plotdata,aes(`Group`,`Uscore`)) +
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic() +
  stat_compare_means(method = "anova",label = "p.format")

write.csv(plotdata,file = "HCC_Subgroup_plotdata_fin_20220521.csv")
write.csv(uptakegene,file = "Uptake_gene.csv")
write.csv(biosyngene,file = "Biosynthesis_gene.csv")