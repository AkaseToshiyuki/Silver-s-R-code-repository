library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(Matrix)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(IOBR)

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

#Load data
rawM <- read.table("Table5A_MOCK.csv",sep = ",",header = F)
rawS <- read.table("Table5B_AntiCD19.csv",sep = ",",header = F)
groupM <- rawM[1:2,]
groupS <- rawS[1:2,]
rawM <- rawM[-1,]
rawS <- rawS[-1,]
rawM <- rawM[which(!duplicated(rawM$V1)),]
rawS <- rawS[which(!duplicated(rawS$V1)),]
rownames(rawM) <- rawM$V1
rownames(rawS) <- rawS$V1
rawM <- rawM[,-1]
rawS <- rawS[,-1]
colnames(rawM) <- rawM[1,]
colnames(rawS) <- rawS[1,]
rawM <- rawM[-1,]
rawS <- rawS[-1,]

rawM1 <- rawM
rawS1 <- rawS

for (i in 1:ncol(rawM)) {
    rawM1[,i] <- as.numeric(rawM[,i])
}
rawM1 <- rawM1[1:nrow(rawM),]
colnames(rawM1) <- colnames(rawM)
rownames(rawM1) <- rownames(rawM)

for (i in 1:ncol(rawS)) {
  rawS1[,i] <- as.numeric(rawS[,i])
}
rawS1 <- rawS1[1:nrow(rawM),]
colnames(rawS1) <- colnames(rawS)
rownames(rawS1) <- rownames(rawS)

rawMtpm <- count2tpm(rawM1,idType = "SYMBOL",org = "hsa",source = "web")
rawStpm <- count2tpm(rawS1,idType = "SYMBOL",org = "hsa",source = "web")

rawMtpm <- as.data.frame(rawMtpm)
rawStpm <- as.data.frame(rawStpm)

save(rawMtpm,file = "Bulkseq_rawMock_group.RData")
save(rawStpm,file = "Bulkseq_rawCD19_group.RData")