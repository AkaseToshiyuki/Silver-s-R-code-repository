library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(clusterProfiler) 
library(ggplot2)


biosyngene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24")
syngene <- c("SREBF2","HMGCR","SQLE")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
CD8gene <- c("CD8A","CD8B")
uptakegene <- c("LDLR", "NPC1", "NPC2")

raw <- read.table("GSE146771_CRC.Leukocyte.Smart-seq2.TPM.txt.gz",sep = " ",header = T)
cell <- read.table("GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz",sep = "\t",header = T)
table(cell$Global_Cluster)
table(cell$Sub_Cluster)
table(cell$Tissue)

Tumorcell <- cell[which(cell$Tissue == "T"),]
cluster <- split(Tumorcell$CellName,Tumorcell$Global_Cluster)

geneused <- uptakegene

B <- raw[geneused,which(colnames(raw) %in% cluster$`B cell`)]
CD4 <- raw[geneused,which(colnames(raw) %in% cluster$`CD4 T cell`)]
CD8 <- raw[geneused,which(colnames(raw) %in% cluster$`CD8 T cell`)]
Epi <- raw[geneused,which(colnames(raw) %in% cluster$"Epithelial cell")]
Fib <- raw[geneused,which(colnames(raw) %in% cluster$Fibroblast)]
ILC <- raw[geneused,which(colnames(raw) %in% cluster$ILC)]
MaC <- raw[geneused,which(colnames(raw) %in% cluster$`Malignant cell`)]
MyC <- raw[geneused,which(colnames(raw) %in% cluster$`Myeloid cell`)]

B["Mean",] <- colMeans(B)
CD4["Mean",] <- colMeans(CD4)
CD8["Mean",] <- colMeans(CD8)
Epi["Mean",] <- colMeans(Epi)
Fib["Mean",] <- colMeans(Fib)
ILC["Mean",] <- colMeans(ILC)
MaC["Mean",] <- colMeans(MaC)
MyC["Mean",] <- colMeans(MyC)

B["Group",] <- "B Cell"
CD4["Group",] <- "CD4 T Cell"
CD8["Group",] <- "CD8 T Cell"
Epi["Group",] <- "Epithelial Cell"
Fib["Group",] <- "Fibroblast"
ILC["Group",] <- "ILC"
MaC["Group",] <- "Malignant Cell"
MyC["Group",] <- "Myeloid Cell"

metadata <- bind_cols(B,CD4,CD8,Epi,Fib,ILC,MaC,MyC)
plotdata <- t(metadata[c("Mean","Group"),]) %>% as.data.frame()
colnames(plotdata)[1] <- "Score"
plotdata$Score <- as.numeric(plotdata$Score)

#write.csv(plotdata,"LXR_Tumor.csv")

ggplot(plotdata,aes(`Group`,`Score`)) + 
  labs(x = "Group", y = "Score", title = "Uptake",theme_classic()) +
  geom_violin(aes(fill = `Score`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20")


