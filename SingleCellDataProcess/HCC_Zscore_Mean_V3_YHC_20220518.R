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

tumorID <- cell[which(cell$Tissue == "Tumor"),]
table(tumorID$Sample)
table(tumorID$celltype_sub)
CD8ID <- tumorID[which(str_sub(tumorID$celltype_sub,1,3) == "CD8"),]
table(CD8ID$celltype_sub)

geneused <- c(acquisitiongene,proliferation)

CD8data <- raw[CD8ID$Barcode,which(colnames(raw) %in% geneused)]
CD8data <- na.omit(CD8data)

scaledata <- as.data.frame(scale(CD8data))

scaledata$PLIscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% proliferation)])

for (i in 1:nrow(scaledata)) {
  if (scaledata$PLIscore[i] >= median(scaledata$PLIscore)) {
    scaledata$Group[i] <- "ProlifHi"
  }else{scaledata$Group[i] <- "ProlifLo"}
}

scaledata$pscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% acquisitiongene)])
scaledata$pscoreM <- sapply(1:nrow(scaledata),function(x){median(unlist(scaledata[x,acquisitiongene]))})

plotdata <- scaledata[,c("pscoreM","pscore","Group","PLIscore")]
plot(plotdata$pscoreM)
summary(plotdata)

ggplot(plotdata,aes(`Group`,`pscore`)) +
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic() +
  stat_compare_means(method = "t.test",label = "p.format")

write.csv(plotdata,file = "HCC_CD8+_plotdata_fin_20220518.csv")
write.csv(proliferation,file = "ProliferationHCC_gene.csv")
write.csv(acquisitiongene,file = "Chol_Acquisition_gene.csv")