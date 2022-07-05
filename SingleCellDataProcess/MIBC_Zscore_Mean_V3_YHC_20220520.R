library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

SREBPgene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
               "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
               "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24","LDLR")
acquisitiongene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                     "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                     "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24","LDLR", "NPC1", "NPC2")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
#LXRgeneN <- c("ABCA1", "ABCG1", "ARL4C", "PLTP", "CETP", "LPL", "APOC1", "APOC2", "APOE", "ABCG5", "ABCG8","NR1H3","NR1H2","MYLIP")
#LXRgeneNP <- c("ABCA1", "ABCG1", "ARL4C", "PLTP", "APOC1", "APOE","MYLIP","NR1H3","NR1H2")
#LXRgeneNPE <- c("ABCA1", "ABCG1", "ARL4C", "PLTP", "CETP", "LPL","APOC1", "APOC2", "APOE", "ABCG5", "ABCG8","NR1H3")
#LXRgeneT <- c("APOE","NR1H2","ABCG1","ABCA1","APOC1")
uptakegene <- c("LDLR", "NPC1", "NPC2")
uptakesubtype <- c("SCARB1","LDLR", "NPC1", "NPC2")
biosyngene <- c("ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24")
CDmarker <- c("CD8A","CD8B","CD4","CD3E")
SREBF2genes <- c("SREBF2","HMGCR","SQLE","LDLR")
functiongene <- c("IFNG", "IL10", "GZMA", "GZMB", "PRF1","TNFSF10", "FASLG", "TNF", "IL2","CD44", "IL2RA")
proliferation <- read.table("ProliferationHCC_gene.csv",header = T,sep = ",")
proliferation <- proliferation[,2]

raw <- read.table("GSE149652_CD8_TIL_droplet_count_matrice.csv.gz",sep = ",",header = T)
save(raw,file = "MIBC_metacount.RData")
#Load data and TPM count
load("MIBC_metacount.RData")
cell <- read.table("GSE149652_CD8_TIL_droplet_cellinfo_matrice.csv.gz",sep = ",",header = T)

cell <- cell[which(cell$tissue == "Tumor" & cell$treatment == "None"),]
cell$index <- str_replace_all(cell$index,"-",".")
rownames(raw) <- raw$index
raw <- raw[,-1]
raw <- raw[,cell$index]

geneinfo <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','transcript_length'),
                  mart = ensembl)
geneinfo <- geneinfo[which(!duplicated(geneinfo$ensembl_gene_id)),]
geneinfo <- geneinfo[which(!duplicated(geneinfo$hgnc_symbol)),]
rownames(geneinfo) <- geneinfo$hgnc_symbol
geneinfo <- geneinfo[rownames(raw),]
geneinfo <- na.omit(geneinfo)
raw <- raw[which(rownames(raw) %in% rownames(geneinfo)),]
geneinfo <- geneinfo[rownames(raw),]

pb <- txtProgressBar(min = 1,max = ncol(raw),style = 3)
for (i in 1:ncol(raw)) {
  raw[,i] <- raw[,i]/geneinfo$transcript_length
  setTxtProgressBar(pb,i)
}

colSums(raw[,1:10])

pb <- txtProgressBar(min = 1,max = ncol(raw),style = 3)
for (i in 1:ncol(raw)) {
  raw[,i] <- raw[,i]/sum(raw[,i])*10^6
  setTxtProgressBar(pb,i)
}

colSums(raw[,1:10])
save(raw,file = "MIBC_TPM_fin.RData")

#Zscore count
load("MIBC_TPM_fin.RData")
geneused <- c(acquisitiongene,proliferation)
raw <- as.data.frame(t(raw))
raw <- raw[which(raw$CD8A > 0 & raw$CD8B > 0 & raw$CD4 == 0),]
raw <- raw[,which(colnames(raw) %in% geneused & colSums(raw) != 0)]
scaledata <- as.data.frame(scale(raw))

scaledata$PLIscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% proliferation)])

for (i in 1:nrow(scaledata)) {
  if (scaledata$PLIscore[i] >= median(scaledata$PLIscore)) {
    scaledata$Group[i] <- "ProlifHi"
  }else{scaledata$Group[i] <- "ProlifLo"}
}

scaledata$pscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% acquisitiongene)])
scaledata$pscoreM <- sapply(1:nrow(scaledata),function(x)
  {median(unlist(scaledata[x,which(colnames(scaledata) %in% acquisitiongene)]))})

plotdata <- scaledata[,c("pscoreM","pscore","Group","PLIscore")]
summary(plotdata)

ggplot(plotdata,aes(`Group`,`pscore`)) +
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  #geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic() +
  stat_compare_means(method = "t.test",label = "p.format")

write.csv(proliferation,file = "ProliferationPDAC_gene.csv")
write.csv(acquisitiongene,file = "Chol_Acquisition_gene.csv")
write.csv(plotdata,file = "PDAC_CD8+_plotdata_fin_20220520.csv")


