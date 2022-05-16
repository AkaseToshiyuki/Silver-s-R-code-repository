library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(DESeq2)
library(dplyr)
library(stringr)

#导入数据并处理
data.count <- read.table("TCGA-BRCA.htseq_counts.tsv.gz",check.names = F,
                         row.names = 1,header = T)
data <- as.matrix(2^data.count-1)
as.character(data[1:100,1:10])
expr <- apply(data, 2, as.integer)
expr[1:4,1:4]
rownames(expr) <- rownames(data)
meta <- as.data.frame(t(expr))
meta[1:4,1:5]

#读入分类信息并给测序数据打tag
clin <- read.csv("BRCA_clinicalMatrix(1).txt",sep = "\t",header = T)
clin [clin == ""] <- NA
pam <- as.data.frame(clin$PAM50Call_RNAseq)
pam[,2] <- clin$sampleID
pam <- dplyr::filter(pam,!is.na(pam$`clin$PAM50Call_RNAseq`))
rownames(pam) <- pam$V2

#筛选具有分类数据的测序结果
same <- meta[which(str_sub(rownames(meta),1,15) %in% rownames(pam)),]
pamt <- pam[which(rownames(pam) %in% str_sub(rownames(same),1,15)),]
same[,60489] <- str_sub(rownames(same),1,15)
same <- same[!duplicated(same$V60489),]
rownames(same) <- same$V60489
same <- same[,-(60484:60489)]
same <- as.matrix(t(same))
same[1:5,1:5]

#DESeq2得到Log2FC
group <- factor(pamt$`clin$PAM50Call_RNAseq`)
dds <- DESeqDataSetFromMatrix(same, DataFrame(group), design= ~ group )
dds <- DESeq(dds)

LumA <- results(dds,contrast = c("group","Basal","LumA"))
LumB <- results(dds,contrast = c("group","Basal","LumB"))
Her2 <- results(dds,contrast = c("group","Basal","Her2"))
Normal <- results(dds,contrast = c("group","Basal","Normal"))

#以上四组仅使用LumA作为演示
dLumA <- as.data.frame(LumA)
dLumA[,"Gene"] <- rownames(dLumA)
dLumA <- dLumA[,c("Gene","log2FoldChange")]

#Gene ID 转换
rownames(dLumA) <- c(1:nrow(dLumA))
dLumA$Gene <- str_sub(dLumA$Gene,1,15)
dLumAgene <- bitr(dLumA[,1],
                  fromType = "ENSEMBL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Hs.eg.db)

#ENTREZID 去重，筛选可以mapping的ENSG ID
dLumAgene <- dLumAgene[!duplicated(dLumAgene$ENTREZID),]
LumA_fin <- dLumA[which(dLumA$Gene  %in% dLumAgene$ENSEMBL),]
LumAgene_fin <- dLumAgene[which(dLumAgene$ENSEMBL %in% LumA_fin$Gene),]
colnames(LumA_fin)[1] <- "ENSEMBL"

#以ENSG ID为桥梁对Log2FC和ENTREZ ID 进行merge
A_fin <- merge(LumA_fin,dLumAgene_fin,by = "ENSEMBL")
A_fin <- A_fin[!duplicated(A_fin$ENTREZID),]
A_fin <- A_fin[,-1]
A_fin <- A_fin[,c("ENTREZID","log2FoldChange")]

#将表达量输出为字符串并命名，最后按照Log2FC从大到小排序
listLumA <- A_fin$log2FoldChange
names(listLumA) <- A_fin$ENTREZID
listLumA <- sort(listLumA,decreasing = T)

#KEGG通路富集并展示结果
keggLumA <- gseKEGG(listLumA,organism = "hsa")
keggLumA
gseaplot(keggLumA, "hsa04658", color.line = "red")
gseaplot(keggLumA, "hsa05140", color.line = "blue")
