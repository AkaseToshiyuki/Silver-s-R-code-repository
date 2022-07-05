library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(ggpubr)

biosyngene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
uptakegene <- c("LDLR", "NPC1", "NPC2")
SREBF2genes <- c("SREBF2","HMGCR","SQLE","LDLR")

cell <- read.table("GSE140228_cell_info_Smartseq2.tsv.gz",sep = "\t",header = T)

load("TPMData_HCC.RData")
TPM[1:5,1:5]
raw <- TPM

tumorID <- cell[which(cell$Tissue == "Tumor"),]
table(tumorID$Sample)
table(tumorID$celltype_sub)
CD8ID <- tumorID[which(str_sub(tumorID$celltype_sub,1,3) == "CD8"),]
table(CD8ID$celltype_sub)

C2 <- CD8ID[which(str_sub(CD8ID$celltype_sub,5,6)== "C2"),]
C4 <- CD8ID[which(str_sub(CD8ID$celltype_sub,5,6)== "C4"),]
C5 <- CD8ID[which(str_sub(CD8ID$celltype_sub,5,6)== "C5"),]
C6 <- CD8ID[which(str_sub(CD8ID$celltype_sub,5,6)== "C6"),]
C7 <- CD8ID[which(str_sub(CD8ID$celltype_sub,5,6)== "C7"),]
C8 <- CD8ID[which(str_sub(CD8ID$celltype_sub,5,6)== "C8"),]
C9 <- CD8ID[which(str_sub(CD8ID$celltype_sub,5,6)== "C9"),]

dtC2 <- t(raw[,which(colnames(raw) %in% C2$Barcode)]) %>% as.data.frame()
dtC4 <- t(raw[,which(colnames(raw) %in% C4$Barcode)]) %>% as.data.frame()
dtC5 <- t(raw[,which(colnames(raw) %in% C5$Barcode)]) %>% as.data.frame()
dtC6 <- t(raw[,which(colnames(raw) %in% C6$Barcode)]) %>% as.data.frame()
dtC7 <- t(raw[,which(colnames(raw) %in% C7$Barcode)]) %>% as.data.frame()
dtC8 <- t(raw[,which(colnames(raw) %in% C8$Barcode)]) %>% as.data.frame()
dtC9 <- t(raw[,which(colnames(raw) %in% C9$Barcode)]) %>% as.data.frame()


dtC2$Group <- "CD8-C2-MKI67"
dtC4$Group <- "CD8-C4-CX3CR1"
dtC5$Group <- "CD8-C5-SELL"
dtC6$Group <- "CD8-C6-GZMK"
dtC7$Group <- "CD8-C7-KLRD1"
dtC8$Group <- "CD8-C8-PDCD1"
dtC9$Group <- "CD8-C9-SLC4A10"

positivegene <- c(biosyngene,uptakegene)
negativegene <- LXRgene

#metadata <- bind_rows(dtC2,dtC4,dtC5,dtC6,dtC7,dtC8,dtC9)

metadata <- bind_rows(dtC2,dtC8)

positivedata <- metadata[,c(positivegene)]
positive <- as.data.frame(sapply(1:(ncol(positivedata)),function(x){positivedata[,x]/mean(positivedata[,x])}))
rownames(positive) <- rownames(positivedata)
colnames(positive) <- colnames(positivedata)

negativedata <- metadata[,c(negativegene)]
negative <- as.data.frame(sapply(1:(ncol(negativedata)),function(x){negativedata[,x]/mean(negativedata[,x])}))
rownames(negative) <- rownames(negativedata)
colnames(negative) <- colnames(negativedata)
negative <- -negative

findata <- cbind(negative,positive)
findata$Score <- rowSums(findata)
findata$Group <- metadata$Group
plotdata <- findata[,c("Group","Score")]

ggplot(plotdata,aes(`Group`,`Score`))+ 
  geom_violin(aes(fill = `Group`))+
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0)+
  geom_jitter(size = 0.2, width = 0.1,shape = 16,alpha = 0.5,color = "grey20")+
  theme_classic()+
  stat_compare_means(method = "t.test",label = "p.format")


plotdata$Ki67 <- metadata$MKI67

ggplot(plotdata,aes(`Score`,`Ki67`))+ 
  geom_point()