library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)


SREBPgene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
               "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
               "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24")
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
proliferation <- read.table("ProliferationBCC_gene.csv",header = T,sep = ",")
proliferation <- proliferation[,2]

raw <- read.table("GSE146771_CRC.Leukocyte.Smart-seq2.TPM.txt.gz",sep = " ",header = T)
raw[1:10,1:10]
raw <- as.data.frame(t(raw))
save(raw,file = "CRC_rawdata_TPM.RData")

#data process
cell <- read.table("GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz",sep = "\t",header = T)
load("CRC_rawdata_TPM.RData")
geneused <- c(uptakegene,biosyngene)

cell <- cell[which(cell$Tissue == "T"),]
raw <- raw[cell$CellName,geneused]

scaledata <- as.data.frame(scale(raw))
scaledata$Group <- cell$Global_Cluster
table(scaledata$Group)

scaledata <- scaledata[which(scaledata$Group !=  "Fibroblast"),]
scaledata <- scaledata[which(scaledata$Group !=  "Epithelial cell"),]
scaledata <- scaledata[which(scaledata$Group !=  "Malignant cell"),]

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

write.csv(plotdata,file = "CRC_Subgroup_plotdata_fin_20220521.csv")
write.csv(uptakegene,file = "Uptake_gene.csv")
write.csv(biosyngene,file = "Biosynthesis_gene.csv")
