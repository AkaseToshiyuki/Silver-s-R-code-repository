library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(Matrix)
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


#load data
barcode.path <- "GSE152938_RAW/GSM4630030_chRCC/barcodes.tsv.gz"
features.path <- "GSE152938_RAW/GSM4630030_chRCC/features.tsv.gz"
matrix.path <- "GSE152938_RAW/GSM4630030_chRCC/matrix.mtx.gz"

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
rawchrcc <- mat
save(rawchrcc,file = "GSE152938_RAW/rawCHRCC.RData")

#count TPM value
load("GSE152938_RAW/rawCHRCC.RData")
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

geneinfo <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','transcript_length'),
                  mart = ensembl)
geneinfo <- geneinfo[which(!duplicated(geneinfo$ensembl_gene_id)),]

rownames(geneinfo) <- geneinfo$ensembl_gene_id
geneinfo <- geneinfo[feature.names$V1,]
geneinfo <- na.omit(geneinfo)
rawchrcc <- rawchrcc[geneinfo$ensembl_gene_id,]
rownames(feature.names) <- feature.names$V1
feature.names <- feature.names[geneinfo$ensembl_gene_id,]

feature.names$Length <- geneinfo$transcript_length
feature.names <- feature.names[which(!duplicated(feature.names$V2)),]
feature.names <- na.omit(feature.names)
rawchrcc <- rawchrcc[feature.names$V1,]
rownames(rawchrcc) <- feature.names$V2
rawchrcc[1:10,1:10]

pb <- txtProgressBar(min = 1,max = ncol(rawchrcc),style = 3)
for (i in 1:ncol(rawchrcc)) {
  rawchrcc[,i] <- rawchrcc[,i]/feature.names$Length
  setTxtProgressBar(pb,i)
}

colSums(rawchrcc[,1:10])

pb <- txtProgressBar(min = 1,max = ncol(rawchrcc),style = 3)
for (i in 1:ncol(rawchrcc)) {
  rawchrcc[,i] <- rawchrcc[,i]/sum(rawchrcc[,i])*10^6
  setTxtProgressBar(pb,i)
}

colSums(rawchrcc[,1:10])
rawchrcc <- as.data.frame(t(rawchrcc))
plot(rawchrcc$CD8A,rawchrcc$CD4)
rawchrcc <- rawchrcc[which(rawchrcc$CD8A > 0 & rawchrcc$CD8B > 0 & rawchrcc$CD4 == 0),]
save(rawchrcc,file = "GSE152938_RAW/TPMchrcc.RData")

#count zscore
load("GSE152938_RAW/TPMchrcc.RData")
load("GSE152938_RAW/TPMCCRCC1.RData")
load("GSE152938_RAW/TPMCCRCC2.RData")
load("GSE152938_RAW/TPMPRCC.RData")
geneused <- c(acquisitiongene,proliferation)

raw <- bind_rows(rawccrcc1,rawccrcc2,rawchrcc,rawprcc)
raw <- raw[,geneused]
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


