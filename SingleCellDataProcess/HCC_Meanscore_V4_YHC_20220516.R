library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(MASS)

biosyngene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24","LDLR", "NPC1", "NPC2")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
uptakegene <- c("LDLR", "NPC1", "NPC2")
SREBF2genes <- c("SREBF2","HMGCR","SQLE","LDLR")
functiongene <- c("IFNG", "IL10", "GZMA", "GZMB", "PRF1","TNFSF10", "FASLG", "TNF", "IL2","CD44", "IL2RA")
proliferation <- read.table("progenefin.txt",header = T,sep = ",")
proliferation <- proliferation[,2]

cell <- read.table("GSE140228_cell_info_Smartseq2.tsv.gz",sep = "\t",header = T)
load("TPMData_HCC.RData")
raw <- t(TPM) %>% as.data.frame()
raw[1:5,1:5]

tumorID <- cell[which(cell$Tissue == "Tumor"),]
table(tumorID$Sample)
table(tumorID$celltype_sub)
CD8ID <- tumorID[which(str_sub(tumorID$celltype_sub,1,3) == "CD8"),]
table(CD8ID$celltype_sub)

geneused <- c(biosyngene,LXRgene)

CD8data <- raw[CD8ID$Barcode,geneused]
CD8data <- na.omit(CD8data)
CD8data[1:10,1:10]
scaledata <- CD8data[,geneused]

scaledata$pscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% biosyngene)])
scaledata$nscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% LXRgene)])
scaledata$score <- scaledata$pscore - scaledata$nscore

plotdata <- scaledata[,c("score","pscore",'nscore')]
plot(plotdata$pscore)
summary(plotdata)


CD8prof <- raw[CD8ID$Barcode,which(colnames(raw) %in% proliferation)]
CD8prof <- na.omit(CD8prof)
progene <- colnames(CD8prof[,which(colnames(CD8prof) %in% proliferation)])

prolifdata <- CD8prof[,progene]
prolifscore <- prolifdata

prolifscore$PLIscore <- rowMeans(prolifscore[,which(colnames(prolifscore) %in% progene)])
plot(prolifscore$PLIscore)

plotdata2 <- cbind(prolifscore$PLIscore,plotdata)
colnames(plotdata2)[1] <- "PLIscore"

ProlifHi <- plotdata2[which(plotdata2$PLIscore > median(plotdata2$PLIscore)),]
ProlifHi$Group <- "ProlifHi"
ProlifLo <- plotdata2[which(plotdata2$PLIscore <= median(plotdata2$PLIscore)),]
ProlifLo$Group <- "ProlifLo"
Prolifbase <- rbind(ProlifHi,ProlifLo)

ggplot(Prolifbase,aes(`Group`,`nscore`)) +
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic() +
  stat_compare_means(method = "t.test",label = "p.format")

