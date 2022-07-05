library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(clusterProfiler) 
library(ggplot2)
library(IOBR)

biosyngene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
uptakegene <- c("LDLR", "NPC1", "NPC2")
SREBF2genes <- c("SREBF2","HMGCR","SQLE")

cell <- read.table("GSE140228_cell_info_Smartseq2.tsv.gz",sep = "\t",header = T)


#从这里开始不要跑！
raw <- read.table("GSE140228_read_counts_Smartseq2.csv.gz",sep = ",",header = T)
colnames(raw) <- str_replace_all(colnames(raw),"\\.","-")

rownames(raw) <- raw[,1]
raw <- raw[,-1]
table(cell$Tissue)
table(cell$celltype_sub)


datatest <- as.matrix(raw)
save(datatest,file = "TPMcount.RData")
#off
load("TPMcount.RData")
TPM <- count2tpm(datatest,idType = "SYMBOL", org = "hsa", source = "web", 
                 effLength = NULL, id = "id",gene_symbol = "symbol", 
                 length = "eff_length", remove_redundancy = "mean")
TPM <- as.data.frame(TPM)
save(TPM,file = "TPMData_HCC.RData")
#从这一行往下开始跑！

load("TPMData_HCC.RData")
TPM[1:5,1:5]
raw <- TPM

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

dttumorCD4 <- t(raw[,which(colnames(raw) %in% tumorCD4$Barcode)]) %>% as.data.frame()
dttumorCD8 <- t(raw[,which(colnames(raw) %in% tumorCD8$Barcode)]) %>% as.data.frame()
dttumorDC <- t(raw[,which(colnames(raw) %in% tumorDC$Barcode)]) %>% as.data.frame()
dttumorB <- t(raw[,which(colnames(raw) %in% tumorB$Barcode)]) %>% as.data.frame()
dttumorMast <- t(raw[,which(colnames(raw) %in% tumorMast$Barcode)]) %>% as.data.frame()
dttumorMono <- t(raw[,which(colnames(raw) %in% tumorMono$Barcode)]) %>% as.data.frame()
dttumorM <- t(raw[,which(colnames(raw) %in% tumorM$Barcode)]) %>% as.data.frame()
dttumorNK <- t(raw[,which(colnames(raw) %in% tumorNK$Barcode)]) %>% as.data.frame()

dttumorCD4$Group <- "CD4"
dttumorCD8$Group <- "CD8"
dttumorDC$Group <- "DC"
dttumorB$Group <- "B"
dttumorMast$Group <- "Mast"
dttumorMono$Group <- "Mono"
dttumorM$Group <- "M"
dttumorNK$Group <- "NK"

usedgene <- biosyngene
metadata <- bind_rows(dttumorCD4,dttumorCD8,dttumorB,dttumorDC,
                      dttumorM,dttumorMast,dttumorMono,dttumorNK)
metadata <- metadata[,c(usedgene,"Group")]
metadata$Score <- rowMeans(metadata[,1:ncol(metadata)-1])
plotdata <- metadata[,(ncol(metadata)-1):ncol(metadata)]


ggplot(plotdata,aes(`Group`,`Score`))+ 
  geom_violin(aes(fill = `Group`))+
  geom_boxplot(width = 0.4,fill = "white",outlier.alpha = 0)+
  geom_jitter(size = 0.2, width = 0.1,shape = 16,alpha = 0.5,color = "grey20")+
  theme_classic()


cd8plotdata <- plotdata[plotdata$Group == "CD8",]
cd8plotdata$PDCD1 <- dttumorCD8$PDCD1

ggplot(cd8plotdata,aes(`Score`,`PDCD1`))+ 
  geom_point()
