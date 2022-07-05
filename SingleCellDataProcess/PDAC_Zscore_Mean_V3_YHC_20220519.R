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


#below is data matrix process
matrix_dir = "GSE155698_RAW/PDAC_TISSUE_9/" 
barcode.path <- paste0(matrix_dir, "filtered_feature_bc_matrix/barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "filtered_feature_bc_matrix/features.tsv.gz")
matrix.path <- paste0(matrix_dir, "filtered_feature_bc_matrix/matrix.mtx.gz")

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
raw9 <- mat
save(raw9,file = "GSE155698_RAW/T9.RData")


load("GSE155698_RAW/T1.RData")
load("GSE155698_RAW/T2.RData")
load("GSE155698_RAW/T3.RData")
load("GSE155698_RAW/T4.RData")
load("GSE155698_RAW/T5.RData")
load("GSE155698_RAW/T6.RData")
load("GSE155698_RAW/T7.RData")
load("GSE155698_RAW/T8.RData")
load("GSE155698_RAW/T9.RData")
rawT1 <- bind_cols(raw1,raw2,raw3,raw4,raw5,raw6,raw7,raw8,raw9)
save(rawT1,file = "GSE155698_RAW/rawT1.RData")

load("GSE155698_RAW/T10.RData")
load("GSE155698_RAW/T11A.RData")
load("GSE155698_RAW/T11B.RData")
load("GSE155698_RAW/T12.RData")
load("GSE155698_RAW/T13.RData")
load("GSE155698_RAW/T15.RData")
load("GSE155698_RAW/T16.RData")
rawT2 <- bind_cols(raw10,raw11A,raw11B,raw12,raw13,raw15,raw16)
save(rawT2,file = "GSE155698_RAW/rawT2.RData")

#TPM count
load("GSE155698_RAW/rawT2.RData")
geneinfo <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','transcript_length'),
                  mart = ensembl)
geneinfo <- geneinfo[which(!duplicated(geneinfo$ensembl_gene_id)),]

rownames(geneinfo) <- geneinfo$ensembl_gene_id
geneinfo <- geneinfo[feature.names$V1,]
geneinfo <- na.omit(geneinfo)
rawT2 <- rawT2[geneinfo$ensembl_gene_id,]
rownames(feature.names) <- feature.names$V1
feature.names <- feature.names[geneinfo$ensembl_gene_id,]

feature.names$Length <- geneinfo$transcript_length
feature.names <- feature.names[which(!duplicated(feature.names$V2)),]
feature.names <- na.omit(feature.names)
rawT2 <- rawT2[feature.names$V1,]
rownames(rawT2) <- feature.names$V2
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

rawT2[1:10,1:10]
rawT2["CD8A",1:10]
colSums(rawT2[1:100])
save(rawT2,file = "GSE155698_RAW/PDAC_TPM2#.RData")
#Get TPM data

load("GSE155698_RAW/PDAC_TPM1#.RData")
rawT1 <- as.data.frame(t(rawT1))
rawT1[1:10,1:10]
plot(rawT1$CD8A,rawT1$CD4)
rawT1 <- rawT1[which(rawT1$CD8A > 0 & rawT1$CD8B > 0 & rawT1$CD4 ==0),]
save(rawT1,file = "GSE155698_RAW/CD8+_TPM_rawT1.RData")

load("GSE155698_RAW/PDAC_TPM2#.RData")
rawT2 <- as.data.frame(t(rawT2))
rawT2[1:10,1:10]
plot(rawT2$CD8A,rawT2$CD4)
rawT2 <- rawT2[which(rawT2$CD8A > 0 & rawT2$CD8B > 0 & rawT2$CD4 ==0),]
save(rawT2,file = "GSE155698_RAW/CD8+_TPM_rawT2.RData")

#Z-score count
geneused <- c(acquisitiongene,proliferation)

raw <- bind_rows(rawT1,rawT2)
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

write.csv(proliferation,file = "ProliferationPDAC_gene.csv")
write.csv(acquisitiongene,file = "Chol_Acquisition_gene.csv")
write.csv(plotdata,file = "PDAC_CD8+_plotdata_fin_20220520.csv")


