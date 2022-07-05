library(dplyr)
library(stringr)
library(Matrix)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
#The upper code just run once is enough


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