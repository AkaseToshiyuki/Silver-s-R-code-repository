library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(Matrix)
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


#below is data matrix process (load 1-18 patients)
matrix_dir = "GSE164690_RAW/" 
barcode.path <- paste0(matrix_dir, "GSM5017070_HN18_CD45p_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "GSM5017070_HN18_CD45p_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "GSM5017070_HN18_CD45p_matrix.mtx.gz")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
mat <- as.matrix(mat) %>% as.data.frame()
raw18 <- mat
save(raw18,file = "GSE164690_RAW/raw18.RData")


rawT1 <- bind_cols(raw1,raw2,raw3,raw4,raw5,raw6)
save(rawT1,file = "GSE164690_RAW/rawT1.RData")
rawT2 <- bind_cols(raw7,raw8,raw9,raw10,raw11,raw12)
save(rawT2,file = "GSE164690_RAW/rawT2.RData")
rawT3 <- bind_cols(raw13,raw14,raw15,raw16,raw17,raw18)
save(rawT3,file = "GSE164690_RAW/rawT3.RData")

#TPM count (rawT1-T3)
load("GSE164690_RAW/rawT3.RData")

geneinfo <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','transcript_length'),
                  mart = ensembl)
geneinfo <- geneinfo[which(!duplicated(geneinfo$ensembl_gene_id)),]

rownames(geneinfo) <- geneinfo$ensembl_gene_id
geneinfo <- geneinfo[feature.names$V1,]
geneinfo <- na.omit(geneinfo)
rawT3 <- rawT3[geneinfo$ensembl_gene_id,]
rownames(feature.names) <- feature.names$V1
feature.names <- feature.names[geneinfo$ensembl_gene_id,]

feature.names$Length <- geneinfo$transcript_length
feature.names <- feature.names[which(!duplicated(feature.names$V2)),]
feature.names <- na.omit(feature.names)
rawT3 <- rawT3[feature.names$V1,]
rownames(rawT3) <- feature.names$V2
rawT3[1:10,1:10]

pb <- txtProgressBar(min = 1,max = ncol(rawT3),style = 3)
for (i in 1:ncol(rawT3)) {
  rawT3[,i] <- rawT3[,i]/feature.names$Length
  setTxtProgressBar(pb,i)
}

colSums(rawT3[,1:10])

pb <- txtProgressBar(min = 1,max = ncol(rawT3),style = 3)
for (i in 1:ncol(rawT3)) {
  rawT3[,i] <- rawT3[,i]/sum(rawT3[,i])*10^6
  setTxtProgressBar(pb,i)
}

colSums(rawT3[,1:10])
rawT3["CD8A",1:10]
save(rawT3,file = "GSE164690_RAW/rawT3_TPM.RData")

load("GSE164690_RAW/rawT1_TPM.RData")
rawT1 <- as.data.frame(t(rawT1))
rawT1[1:10,1:10]
plot(rawT1$CD8A,rawT1$CD4)
rawT1 <- rawT1[which(rawT1$CD8A > 0 & rawT1$CD8B > 0 & rawT1$CD4 == 0),]
save(rawT1,file = "GSE164690_RAW/T1_TPM_fin.RData")

load("GSE164690_RAW/rawT2_TPM.RData")
rawT2 <- as.data.frame(t(rawT2))
rawT2[1:10,1:10]
plot(rawT2$CD8A,rawT2$CD4)
rawT2 <- rawT2[which(rawT2$CD8A > 0 & rawT2$CD8B > 0 & rawT2$CD4 == 0),]
save(rawT2,file = "GSE164690_RAW/T2_TPM_fin.RData")

load("GSE164690_RAW/rawT3_TPM.RData")
rawT3 <- as.data.frame(t(rawT3))
rawT3[1:10,1:10]
plot(rawT3$CD8A,rawT3$CD4)
rawT3 <- rawT3[which(rawT3$CD8A > 0 & rawT3$CD8B > 0 & rawT3$CD4 == 0),]
save(rawT3,file = "GSE164690_RAW/T3_TPM_fin.RData")

#Zscore count
load("GSE164690_RAW/T1_TPM_fin.RData")
load("GSE164690_RAW/T2_TPM_fin.RData")
load("GSE164690_RAW/T3_TPM_fin.RData")
geneused <- c(acquisitiongene,proliferation)

raw <- bind_rows(rawT1,rawT2,rawT3)
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

write.csv(proliferation,file = "ProliferationNHSCC_gene.csv")
write.csv(acquisitiongene,file = "Chol_Acquisition_gene.csv")
write.csv(plotdata,file = "NHSCC_CD8+_plotdata_fin_20220520.csv")


