library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)

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

#read data cycle
raw12 <- read.table("GSE183904_RAW/GSM5573477_sample12.csv",sep = ",",header = T)
rownames(raw12) <- raw12$X
raw12 <- raw12[,-1]

raw13 <- read.table("GSE183904_RAW/GSM5573478_sample13.csv",sep = ",",header = T)
rownames(raw13) <- raw13$X
raw13 <- raw13[,-1]

raw14 <- read.table("GSE183904_RAW/GSM5573479_sample14.csv",sep = ",",header = T)
rownames(raw14) <- raw14$X
raw14 <- raw14[,-1]

raw15 <- read.table("GSE183904_RAW/GSM5573480_sample15.csv",sep = ",",header = T)
rownames(raw15) <- raw15$X
raw15 <- raw15[,-1]

raw16 <- read.table("GSE183904_RAW/GSM5573481_sample16.csv",sep = ",",header = T)
rownames(raw16) <- raw16$X
raw16 <- raw16[,-1]

raw17 <- read.table("GSE183904_RAW/GSM5573482_sample17.csv",sep = ",",header = T)
rownames(raw17) <- raw17$X
raw17 <- raw17[,-1]

raw18 <- read.table("GSE183904_RAW/GSM5573483_sample18.csv",sep = ",",header = T)
rownames(raw18) <- raw18$X
raw18 <- raw18[,-1]

rawT2 <- bind_cols(raw12,raw13,raw14,raw15,raw16,raw17,raw18)
save(rawT2,file = "GSE183904_RAW/rawT2.RData")

rm(list = ls())
gc()

load("GSE183904_RAW/rawT2.RData")
#TPM count
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

feature.names <- as.data.frame(rownames(rawT2))
colnames(feature.names)[1] <- "SYMBOL"

geneinfo <- getBM(attributes=c('hgnc_symbol','transcript_length'),
                  mart = ensembl)
geneinfo <- geneinfo[which(!duplicated(geneinfo$hgnc_symbol)),]

rownames(geneinfo) <- geneinfo$hgnc_symbol
geneinfo <- geneinfo[feature.names$SYMBOL,]
geneinfo <- na.omit(geneinfo)

rawT2 <- rawT2[geneinfo$hgnc_symbol,]
rownames(feature.names) <- feature.names$SYMBOL
feature.names <- as.data.frame(feature.names[geneinfo$hgnc_symbol,])

feature.names$Length <- geneinfo$transcript_length
feature.names <- feature.names[which(!duplicated(feature.names$`feature.names[geneinfo$hgnc_symbol, ]`)),]
feature.names <- na.omit(feature.names)
rawT2 <- rawT2[feature.names$`feature.names[geneinfo$hgnc_symbol, ]`,]
rawT2[1:10,1:10]

pb <- txtProgressBar(min = 1,max = ncol(rawT2),style = 3)
for (i in 1:ncol(rawT2)) {
  rawT2[,i] <- rawT2[,i]/feature.names$Length
  setTxtProgressBar(pb,i)
}

colSums(rawT2[,1:10])

pb <- txtProgressBar(min = 1,max = ncol(rawT2),style = 3)
for (i in 1:ncol(rawT2)) {
  rawT2[,i] <- rawT2[,i]/sum(rawT2[,i])*10^6
  setTxtProgressBar(pb,i)
}

colSums(rawT2[,1:10])
rawT2["CD8A",1:10]
save(rawT2,file = "GSE183904_RAW/rawT2_TPM.RData")

#data filter set and zscore count
geneused <- c(acquisitiongene,proliferation,"CD8A","CD8B","CD4")

load("GSE183904_RAW/rawT1_TPM.RData")
rawT1 <- rawT1[geneused,]
rawT1 <- as.data.frame(t(rawT1))
plot(rawT1$CD8A,rawT1$CD4)
rawT1 <- rawT1[which(rawT1$CD8A > 0 & rawT1$CD8B > 0 & rawT1$CD4 ==0),]

load("GSE183904_RAW/rawT2_TPM.RData")
rawT2 <- rawT2[geneused,]
rawT2 <- as.data.frame(t(rawT2))
plot(rawT2$CD8A,rawT2$CD4)
rawT2 <- rawT2[which(rawT2$CD8A > 0 & rawT2$CD8B > 0 & rawT2$CD4 ==0),]

load("GSE183904_RAW/rawT3_TPM.RData")
rawT3 <- rawT3[geneused,]
rawT3 <- as.data.frame(t(rawT3))
plot(rawT3$CD8A,rawT3$CD4)
rawT3 <- rawT3[which(rawT3$CD8A > 0 & rawT3$CD8B > 0 & rawT3$CD4 ==0),]

load("GSE183904_RAW/rawT4_TPM.RData")
rawT4 <- rawT4[geneused,]
rawT4 <- as.data.frame(t(rawT4))
plot(rawT4$CD8A,rawT4$CD4)
rawT4 <- rawT4[which(rawT4$CD8A > 0 & rawT4$CD8B > 0 & rawT4$CD4 ==0),]

load("GSE183904_RAW/rawT5_TPM.RData")
rawT5 <- rawT5[geneused,]
rawT5 <- as.data.frame(t(rawT5))
plot(rawT5$CD8A,rawT5$CD4)
rawT5 <- rawT5[which(rawT5$CD8A > 0 & rawT5$CD8B > 0 & rawT5$CD4 ==0),]

raw <- bind_rows(rawT1,rawT2,rawT3,rawT4,rawT5)
raw <- raw[,which(colnames(raw) %in% geneused)]
raw <- raw[,which(colSums(raw) != 0)]

scaledata <- as.data.frame(scale(raw))

scaledata$PLIscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% proliferation)])

for (i in 1:nrow(scaledata)) {
  if (scaledata$PLIscore[i] >= median(scaledata$PLIscore)) {
    scaledata$Group[i] <- "ProlifHi"
  }else{scaledata$Group[i] <- "ProlifLo"}
}

scaledata$pscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% acquisitiongene)])
scaledata$pscoreM <- sapply(1:nrow(scaledata),function(x){median(unlist(scaledata[x,acquisitiongene]))})

plotdata <- scaledata[,c("pscoreM","pscore","Group","PLIscore")]
summary(plotdata)

ggplot(plotdata,aes(`Group`,`pscore`)) +
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  #geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic() +
  stat_compare_means(method = "t.test",label = "p.format")

write.csv(proliferation,file = "ProliferationGC_gene.csv")
write.csv(acquisitiongene,file = "Chol_Acquisition_gene.csv")
write.csv(plotdata,file = "GC_CD8+_plotdata_fin_20220521.csv")
