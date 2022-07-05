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
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
uptakegene <- c("LDLR", "NPC1", "NPC2")
SREBF2genes <- c("SREBF2","HMGCR","SQLE")
CD8gene <- c("CD8A","CD8B")

raw <- read.table("GSE146771_CRC.Leukocyte.Smart-seq2.TPM.txt.gz",sep = " ",header = T)
cell <- read.table("GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz",sep = "\t",header = T)
table(cell$Global_Cluster)
table(cell$Sub_Cluster)
table(cell$Tissue)
cell <- cell[which(cell$Tissue == "T" & str_sub(cell$Sub_Cluster,6,8) == "CD8"),]
table(cell$Sub_Cluster)

cluster <- split(cell$CellName,cell$Sub_Cluster)

geneused <- uptakegene

ht12 <- raw[geneused,which(colnames(raw) %in% cluster$`hT12_CD8-LEF1`)]
ht13 <- raw[geneused,which(colnames(raw) %in% cluster$`hT13_CD8-GPR183`)]
ht14 <- raw[geneused,which(colnames(raw) %in% cluster$`hT14_CD8-CX3CR1`)]
ht15 <- raw[geneused,which(colnames(raw) %in% cluster$`hT15_CD8-GZMK`)]
ht16 <- raw[geneused,which(colnames(raw) %in% cluster$`hT16_CD8-CD6`)]
ht17 <- raw[geneused,which(colnames(raw) %in% cluster$`hT17_CD8-CD160`)]
ht18 <- raw[geneused,which(colnames(raw) %in% cluster$`hT18_CD8-LAYN`)]

ht12["Mean",] <- colMeans(ht12)
ht13["Mean",] <- colMeans(ht13)
ht14["Mean",] <- colMeans(ht14)
ht15["Mean",] <- colMeans(ht15)
ht16["Mean",] <- colMeans(ht16)
ht17["Mean",] <- colMeans(ht17)
ht18["Mean",] <- colMeans(ht18)

ht12["Group",] <- "CD8-LEF1"
ht13["Group",] <- "CD8-GPR183"
ht14["Group",] <- "CD8-CX3CR1"
ht15["Group",] <- "CD8-GZMK"
ht16["Group",] <- "CD8-CD6"
ht17["Group",] <- "CD8-CD160"
ht18["Group",] <- "CD8-LAYN"


metadata <- bind_cols(ht12,ht13,ht14,ht15,ht16,ht17,ht18)
plotdata <- t(metadata[c("Mean","Group"),]) %>% as.data.frame()
colnames(plotdata)[1] <- "Score"
plotdata$Score <- as.numeric(plotdata$Score)

write.csv(plotdata,"Biosynthesis_CD8group.csv")

ggplot(plotdata,aes(`Group`,`Score`)) + 
  labs(x = "Group", y = "Score", title = "LXR Pathway",theme_classic()) +
  geom_violin(aes(fill = `Score`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20")