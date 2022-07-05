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
proliferation <- read.table("progenefin.txt",header = T,sep = ",")
proliferation <- proliferation[,2]

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
#upper just run once

#From here is data treat
matrix_dir = "GSE151511/ac24/" 
barcode.path <- paste0(matrix_dir, "ac24_barcodes.tsv")
features.path <- paste0(matrix_dir, "ac24_genes.tsv")
matrix.path <- paste0(matrix_dir, "ac24_matrix.mtx")

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
raw <- mat
rm(mat)
raw[1:5,1:5]

geneinfo <- getBM(attributes=c('ensembl_gene_id','transcript_length'),
                  mart = ensembl)
geneinfo <- geneinfo[which(!duplicated(geneinfo$ensembl_gene_id)),]
rownames(geneinfo) <- geneinfo$ensembl_gene_id
geneinfo <- geneinfo[feature.names$V1,]
geneinfo <- na.omit(geneinfo)
raw <- raw[geneinfo$ensembl_gene_id,]
rownames(feature.names) <- feature.names$V1
feature.names <- feature.names[geneinfo$ensembl_gene_id,]

feature.names$Length <- geneinfo$transcript_length
feature.names <- feature.names[which(!duplicated(feature.names$V2)),]
feature.names <- na.omit(feature.names)
raw <- raw[feature.names$V1,]
rownames(raw) <- feature.names$V2

pb <- txtProgressBar(min = 1,max = ncol(raw),style = 3)
for (i in 1:ncol(raw)) {
  raw[,i] <- raw[,i]/feature.names$Length
  setTxtProgressBar(pb,i)
}

colSums(raw[,1:10])

pb <- txtProgressBar(min = 1,max = ncol(raw),style = 3)
for (i in 1:ncol(raw)) {
  raw[,i] <- raw[,i]/sum(raw[,i])*10^6
  setTxtProgressBar(pb,i)
}

raw[1:10,1:10]
raw["CD8A",1:10]
sum(raw[,1])
raw_ac24 <- as.data.frame(t(raw))
raw_ac24$Sample <- "ac24"
save(raw_ac24,file = "CART_NM_TPM_20220516/ac24_TPM.RData")
rm(raw_ac24)
rm(raw)
gc()

load("rawf24.RData")

rawf24 <- bind_rows(raw_ac21,raw_ac22,raw_ac23,raw_ac24)
save(rawf24,file = "CART_NM_TPM_20220516/rawf24.RData")

geneused <- c(biosyngene,LXRgeneN,proliferation,functiongene,CDmarker,"Sample")
load("CART_NM_TPM_20220516/rawf24.RData")
data24 <- rawf24[,which(colnames(rawf24) %in% geneused)]
save(data24,file = "CART_NM_TPM_20220516/data24.RData")

geneused <- c(biosyngene,LXRgeneN,proliferation,functiongene,CDmarker,"Sample")
load("CART_NM_TPM_20220516/rawf6.RData")
data6 <- rawf6[,which(colnames(rawf6) %in% geneused)]
save(data6,file = "CART_NM_TPM_20220516/data6.RData")

load("CART_NM_TPM_20220516/data6.RData")
listdata <- bind_rows(data6,data12,data16,data20,data24)
table(listdata$Sample)
save(listdata,file = "CART_Listdata.RData")
#upper is data treat
#Please adjust the code upper to fit your work environment
#DO NOT JUST USE THEM!

#Universal treat from HERE!
#Score acquire from here
load("CART_Listdata.RData")
CD3Esum <- summary(listdata$CD3E)
CD8Asum <- summary(listdata$CD8A)
CD8Bsum <- summary(listdata$CD8B)
CD4sum <- summary(listdata$CD4)

listdata <- listdata[which(listdata$Sample != "ac06"),]
listdata <- listdata[which(listdata$CD3E > 0),]

scaledata <- listdata[,which(colSums(listdata[,1:(ncol(listdata)-1)]) != 0)]
scaledata[1:10,1:10]
scaledata <- as.data.frame(scale(scaledata))
scaledata$Sample <- listdata$Sample

scaledata[1:5,SREBPgene]
scaledata[1:5,LXRgene]
table(scaledata$Sample)

CRlist <- c("ac01","ac05","ac07","ac08","ac09","ac10","ac12","ac14","ac16")
PDlist <- c("ac02","ac03","ac04","ac11","ac13","ac15","ac17","ac18","ac19","ac20","ac21","ac22","ac23","ac24")

for (i in 1:nrow(scaledata)) {
  if (scaledata$Sample[i] %in% CRlist) {
    scaledata$Group[i] <- "CR"}
  else{
    scaledata$Group[i] <- "PD"}  
}

scaledata$pscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% SREBPgene)])
scaledata$nscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% LXRgene)])
scaledata$score <- scaledata$pscore - scaledata$nscore

scaledata$pscoreM <- sapply(1:nrow(scaledata),function(x){median(unlist(scaledata[x,SREBPgene]))})
scaledata$nscoreM <- sapply(1:nrow(scaledata),function(x){median(unlist(scaledata[x,LXRgene]))})
scaledata$scoreM <- scaledata$pscoreM - scaledata$nscoreM

"
tmp <- sapply(1:nrow(scaledata),function(x){which(scaledata[x,] == scaledata$nscoreM[x], arr.ind = TRUE)})
tmp <- as.data.frame(t(tmp))
rownames(tmp) <- rownames(scaledata)
tmp$gene <- sapply(1:nrow(tmp),function(x){colnames(scaledata[tmp$V3[x]])})
table(tmp$gene)
"

plotdata <- scaledata[,c("pscoreM","nscoreM","scoreM","score","pscore",'nscore',"Sample","Group")]
summary(plotdata)
#write.csv(plotdata,file = "CART_NM_plotdata.csv")
#plotdata2 <- plotdata[which(plotdata$nscoreM > 0.3),]

#write.csv(plotdata,file = "plotdt_fin_20220517_CD3E_scale.csv")
ggplot(plotdata,aes(`Group`,`score`)) +
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  #geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic() +
  stat_compare_means(method = "t.test",label = "p.format")

#write.csv(plotdata,"CART_plotdata_fin_20220518.csv")

plot(plotdata$nscore)
plotdata1 <- plotdata[which(plotdata$score >= -5),]
mean(plotdata1$score[which(plotdata1$Group == "CR")])
