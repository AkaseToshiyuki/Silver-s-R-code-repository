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
proliferation <- read.table("ProliferationBCC_gene.csv",header = T,sep = ",")
proliferation <- proliferation[,2]

#data save
raw <- read.table("bcc_scRNA_counts.txt",sep = "\t",header = T)
save(raw,file = "raw_BCC.RData")
raw[1:5,1:5]

#data treat
load("raw_BCC.RData")
raw[1:5,1:5]


searchAttributes(mart = ensembl,pattern = "hgnc")
geneinfo <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','transcript_length'),mart = ensembl)
geneinfo <- geneinfo[which(!duplicated(geneinfo$ensembl_gene_id)),]
geneinfo <- geneinfo[which(!duplicated(geneinfo$hgnc_symbol)),]
rownames(geneinfo) <- geneinfo$hgnc_symbol
geneinfoused <- geneinfo[which(rownames(geneinfo) %in% rownames(raw)),]
raw <- raw[geneinfoused$hgnc_symbol,]
#Save the work frame
save(geneinfoused,file = "BCC_RAWdata_geneinfo.RData")
save(raw,file = "BCC_filtered_RAWdata.RData")

#Load the work frame
load("BCC_RAWdata_geneinfo.RData")
load("BCC_filtered_RAWdata.RData")

pb <- txtProgressBar(min = 1,max = ncol(raw),style = 3)
for (i in 1:ncol(raw)) {
  raw[,i] <- raw[,i]/geneinfoused$transcript_length
  setTxtProgressBar(pb,i)
}

colSums(raw[,1:10])

pb <- txtProgressBar(min = 1,max = ncol(raw),style = 3)
for (i in 1:ncol(raw)) {
  raw[,i] <- raw[,i]/sum(raw[,i])*10^6
  setTxtProgressBar(pb,i)
}

colSums(raw[,1:10])
raw[1:10,1:10]
save(raw,file = "BCC_TPMData.RData")
#Save the TPM data

load("BCC_TPMData.RData")
raw[1:10,1:10]
geneused <- c(biosyngene,functiongene,LXRgene,proliferation)
raw <- raw[geneused,]
raw <- as.data.frame(t(raw))
save(raw,file = "BCC_GenesTPM.RData")

#CD8 subtype analysis
load("BCC_GenesTPM.RData")
cell <- read.table("GSE123813_bcc_all_metadata.txt.gz",sep = "\t",header = T)
table(cell[,3:5])
cell <- cell[which(str_sub(cell$cluster,1,3) == "CD8" & 
                     cell$treatment == "pre" & cell$sort == "CD45+ CD3+"),]
table(cell[,3:5])

CD8data <- raw[which(rownames(raw) %in% cell$cell.id),]
CD8data <- CD8data[,which(colSums(CD8data) != 0)]
CD8data[1:10,acquisitiongene]

proliferationBCC <- colnames(CD8data[,which(colnames(CD8data) %in% proliferation)])
CD8data[1:10,proliferationBCC]
scaledata <- as.data.frame(scale(CD8data))

scaledata$PLIscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% proliferationBCC)])


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

write.csv(proliferationBCC,file = "ProliferationBCC_gene.csv")
write.csv(acquisitiongene,file = "Chol_Acquisition_gene.csv")
write.csv(plotdata,file = "BCC_CD8+_plotdata_fin_20220518.csv")
