library(dplyr)
library(stringr)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)

#raw <- read.table("GSE125881_raw.expMatrix.csv.gz",header = T,sep = ",")
#save(raw,file = "GSE125881_CART_raw.RData")
load("GSE125881_CART_raw.RData")
raw[1:5,1:5]

rownames(raw) <- raw[,1]
raw <- raw[,-1]

rawgene <- rownames(raw) %>% as.data.frame()
gene.df <- bitr(rawgene[,1], fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),OrgDb = org.Hs.eg.db)
head(gene.df)

raw <- raw[gene.df$SYMBOL,]
raw <- cbind(gene.df,raw)
raw <- raw[which(!duplicated(raw$ENSEMBL)),] 
rownames(raw) <- raw$ENSEMBL
gene.df <- raw[,c("SYMBOL","ENSEMBL","ENTREZID")]
raw <- raw[,-(1:3)]

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneinfo <- getBM(attributes=c('ensembl_gene_id','transcript_length'),
                  mart = ensembl)
geneinfo <- geneinfo[which(!duplicated(geneinfo$ensembl_gene_id)),]
rownames(geneinfo) <- geneinfo$ensembl_gene_id
geneinfo <- geneinfo[rownames(raw),]
gene.df$Length <- geneinfo$transcript_length
gene.df <- gene.df[which(!duplicated(gene.df$SYMBOL)),]
gene.df <- na.omit(gene.df)
raw <- raw[gene.df$ENSEMBL,]

pb <- txtProgressBar(min = 1,max = ncol(raw),style = 3)
for (i in 1:ncol(raw)) {
  raw[,i] <- raw[,i]/gene.df$Length
  setTxtProgressBar(pb,i)
}

colSums(raw[,1:10])

pb <- txtProgressBar(min = 1,max = ncol(raw),style = 3)
for (i in 1:ncol(raw)) {
  raw[,i] <- raw[,i]/sum(raw[,i])*10^6
  setTxtProgressBar(pb,i)
}

raw[1:10,1:10]
rownames(raw) <- gene.df$SYMBOL
save(raw,file = "CART_TPMvalue.RData")

