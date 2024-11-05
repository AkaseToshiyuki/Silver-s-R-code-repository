library(DESeq2)
library(biomaRt)
library(stringr)
library(tibble)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpmisc)
library(msigdbr)
library(org.Hs.eg.db)

#TPM
E28z <- read.csv("Others/XXY_CART/20230404_NFAT/E_readsCount(1).csv")
N28Z <- read.csv("Others/XXY_CART/20230404_NFAT/z28_readsCount(1).csv")
B6I <- read.csv("Others/XXY_CART/20230404_NFAT/B_readsCount(1).csv")

counts <- bind_cols(N28Z,E28z[,2:16],B6I[,2:16])
counts$X <- str_extract(counts$X,pattern = ".*(?=\\.)")
counts <- counts[which(!duplicated(counts$X)),]
rownames(counts) <- counts$X
counts <- counts[,2:46]

load("GeneinfoENSG.RData")

geneinfoTPMid <- geneinfo[which(duplicated(geneinfo$ensembl_gene_id)),"ensembl_gene_id"]
geneinfoTPM <- geneinfo[which(geneinfo$ensembl_gene_id %in% geneinfoTPMid),]
geneinfoTPM2 <- geneinfo[which(!(geneinfo$ensembl_gene_id %in% geneinfoTPM$ensembl_gene_id)),]
genenamedup <- geneinfoTPM[which(!duplicated(geneinfoTPM$hgnc_symbol)),"ensembl_gene_id"]

pb <- txtProgressBar(min = 1,max = length(genenamedup),style = 3)
for (i in 1:length(genenamedup)) {
  geneinfoTPM[which(geneinfoTPM$ensembl_gene_id == genenamedup[i]),"transcript_length"] <- 
    max(geneinfoTPM[which(geneinfoTPM$ensembl_gene_id == genenamedup[i]),"transcript_length"])
  setTxtProgressBar(pb,i)
}
geneinfoTPM <- geneinfoTPM[which(!duplicated(geneinfoTPM$ensembl_gene_id)),]
geneinfofin <- bind_rows(geneinfoTPM,geneinfoTPM2)

counts <- add_column(counts,.before = "z28ctrl1",ensembl_gene_id = rownames(counts))
data <- merge(counts,geneinfofin)

pb <- txtProgressBar(min = 2,max = ncol(data)-2,style = 3)
for (i in 2:(ncol(data)-2)) {
  data[,i] <- data[,i]/data$transcript_length
  setTxtProgressBar(pb,i)
}

colSums(data[,2:(ncol(data)-2)])

pb <- txtProgressBar(min = 1,max = (ncol(data)-2),style = 3)
for (i in 2:(ncol(data)-2)) {
  data[,i] <- data[,i]/sum(data[,i])*10^6
  setTxtProgressBar(pb,i)
}

colSums(data[,2:(ncol(data)-2)])
data <- data[which(!duplicated(data$hgnc_symbol)),]
rownames(data) <- data$hgnc_symbol
data <- data[,2:46]

SLE <- read.csv("Speed_Limited_enzyme.csv",header = F)
SLE <- SLE$V1
plotdata <- data[SLE,]
plotdata <- na.omit(plotdata)
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_Speedlimited_Enzyme_Expression.pdf",width = 10,height = 6)
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,ncol(plotdata))]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_boxplot(aes(fill = Groupfin))+
          labs(title = coln[i]))
}
dev.off()

#Glycolysis Glycogen TCAcycle
E28z <- read.csv("Others/XXY_CART/20230404_NFAT/E_readsCount(1).csv")
N28Z <- read.csv("Others/XXY_CART/20230404_NFAT/z28_readsCount(1).csv")
B6I <- read.csv("Others/XXY_CART/20230404_NFAT/B_readsCount(1).csv")

load("GeneinfoENSG.RData")
geneinfo <- geneinfo[which(!duplicated(geneinfo$hgnc_symbol)),]
humanhallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, ensembl_gene)

counts <- bind_cols(N28Z,E28z[,2:16],B6I[,2:16])
counts$X <- str_extract(counts$X,pattern = ".*(?=\\.)")
counts <- counts[which(!duplicated(counts$X)),]
rownames(counts) <- counts$X
counts <- counts[,2:46]

coldata <- as.data.frame(colnames(counts))
coldata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
coldata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
coldata$Groupfin <- paste(coldata$Group,coldata$Round,sep = "_")
colnames(coldata)[1] <- "names"
rownames(coldata) <- coldata$names

counts <- round(counts)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Groupfin)
dds <- DESeq(dds)
Round <- c("Ctrl","R1","R2","R3","R5")

#Glycolysis
gmtgene <- bind_rows(read.gmt("GMTs/HALLMARK_GLYCOLYSIS.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/KEGG_GLYCOLYSIS_GLUCONEOGENESIS.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/REACTOME_GLYCOLYSIS.v2022.1.Hs.gmt"))

pdf(file = "XXY_20230910_Glycolysis.pdf",height = 10,width = 14)
for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs E28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs 28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("E28Z",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()

#TCACycle
gmtgene <- bind_rows(read.gmt("GMTs/KEGG_CITRATE_CYCLE_TCA_CYCLE.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/REACTOME_PYRUVATE_METABOLISM_AND_CITRIC_ACID_TCA_CYCLE.v2023.1.Hs.gmt"))

pdf(file = "XXY_20230910_TCA_Cycle.pdf",height = 10,width = 14)
for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs E28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs 28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("E28Z",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()

#Glycogen
gmtgene <- bind_rows(read.gmt("GMTs/GOBP_POSITIVE_REGULATION_OF_GLYCOGEN_METABOLIC_PROCESS.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/REACTOME_GLYCOGEN_METABOLISM.v2023.1.Hs.gmt"))

pdf(file = "XXY_20230910_Glycogen_Metabolism.pdf",height = 10,width = 14)
for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs E28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs 28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("E28Z",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()

#Merged Glycolysis and TCA cycle
gmtgene <- bind_rows(read.gmt("GMTs/KEGG_CITRATE_CYCLE_TCA_CYCLE.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/KEGG_GLYCOLYSIS_GLUCONEOGENESIS.v2023.1.Hs.gmt"))
gmtgene$term <- "KEGG_GLYCOLYSIS_AND_TCA_CYCLE"

pdf(file = "XXY_20230911_Glycolysis_and_TCA.pdf",height = 10,width = 14)
for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs E28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs 28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("E28Z",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()

#TCA Single Gene
E28z <- read.csv("Others/XXY_CART/20230404_NFAT/E_readsCount(1).csv")
N28Z <- read.csv("Others/XXY_CART/20230404_NFAT/z28_readsCount(1).csv")
B6I <- read.csv("Others/XXY_CART/20230404_NFAT/B_readsCount(1).csv")

counts <- bind_cols(N28Z,E28z[,2:16],B6I[,2:16])
counts$X <- str_extract(counts$X,pattern = ".*(?=\\.)")
counts <- counts[which(!duplicated(counts$X)),]
rownames(counts) <- counts$X
counts <- counts[,2:46]

load("GeneinfoENSG.RData")

geneinfoTPMid <- geneinfo[which(duplicated(geneinfo$ensembl_gene_id)),"ensembl_gene_id"]
geneinfoTPM <- geneinfo[which(geneinfo$ensembl_gene_id %in% geneinfoTPMid),]
geneinfoTPM2 <- geneinfo[which(!(geneinfo$ensembl_gene_id %in% geneinfoTPM$ensembl_gene_id)),]
genenamedup <- geneinfoTPM[which(!duplicated(geneinfoTPM$hgnc_symbol)),"ensembl_gene_id"]

pb <- txtProgressBar(min = 1,max = length(genenamedup),style = 3)
for (i in 1:length(genenamedup)) {
  geneinfoTPM[which(geneinfoTPM$ensembl_gene_id == genenamedup[i]),"transcript_length"] <- 
    max(geneinfoTPM[which(geneinfoTPM$ensembl_gene_id == genenamedup[i]),"transcript_length"])
  setTxtProgressBar(pb,i)
}
geneinfoTPM <- geneinfoTPM[which(!duplicated(geneinfoTPM$ensembl_gene_id)),]
geneinfofin <- bind_rows(geneinfoTPM,geneinfoTPM2)

counts <- add_column(counts,.before = "z28ctrl1",ensembl_gene_id = rownames(counts))
data <- merge(counts,geneinfofin)

pb <- txtProgressBar(min = 2,max = ncol(data)-2,style = 3)
for (i in 2:(ncol(data)-2)) {
  data[,i] <- data[,i]/data$transcript_length
  setTxtProgressBar(pb,i)
}

colSums(data[,2:(ncol(data)-2)])

pb <- txtProgressBar(min = 1,max = (ncol(data)-2),style = 3)
for (i in 2:(ncol(data)-2)) {
  data[,i] <- data[,i]/sum(data[,i])*10^6
  setTxtProgressBar(pb,i)
}

colSums(data[,2:(ncol(data)-2)])
data <- data[which(!duplicated(data$hgnc_symbol)),]
rownames(data) <- data$hgnc_symbol
data <- data[,2:46]

TCA <- read.gmt("GMTs/KEGG_CITRATE_CYCLE_TCA_CYCLE.v2023.1.Hs.gmt")
TCA <- TCA$gene
plotdata <- data[TCA,]
plotdata <- na.omit(plotdata)
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_20230911_TCA_Single_gene.pdf",width = 10,height = 6)
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,ncol(plotdata))]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_boxplot(aes(fill = Groupfin))+
          labs(title = coln[i]))
}
dev.off()

gly <- read.gmt("GMTs/KEGG_GLYCOLYSIS_GLUCONEOGENESIS.v2023.1.Hs.gmt")
gly <- gly$gene
plotdata <- data[gly,]
plotdata <- na.omit(plotdata)
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_20230911_Glycolysis_Single_gene.pdf",width = 10,height = 6)
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,ncol(plotdata))]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_boxplot(aes(fill = Groupfin))+
          labs(title = coln[i]))
}
dev.off()

#Oxphos
gmtgene <- read.gmt("GMTs/KEGG_OXIDATIVE_PHOSPHORYLATION.v2023.1.Hs.gmt")

pdf(file = "XXY_20230911_OXPHOS.pdf",height = 10,width = 14)
for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs E28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs 28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("E28Z",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()

OXP <- read.gmt("GMTs/KEGG_OXIDATIVE_PHOSPHORYLATION.v2023.1.Hs.gmt")
OXP <- OXP$gene
plotdata <- data[OXP,]
plotdata <- na.omit(plotdata)
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_20230911_OXPHOS_Single_gene.pdf",width = 10,height = 6)
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,ncol(plotdata))]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_boxplot(aes(fill = Groupfin))+
          labs(title = coln[i]))
}
dev.off()

FAO <- read.gmt("GMTs/FATTY_ACID_OXIDATION.v2023.1.Hs.gmt")
FAO <- FAO$gene
plotdata <- data[FAO,]
plotdata <- na.omit(plotdata)
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_20230913_FAO_Single_gene.pdf",width = 10,height = 6)
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,ncol(plotdata))]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_boxplot(aes(fill = Groupfin))+
          labs(title = coln[i]))
}
dev.off()

GLU <- read.gmt("GMTs/GOBP_GLUTAMINE_FAMILY_AMINO_ACID_CATABOLIC_PROCESS.v2023.1.Hs.gmt")
GLU <- GLU$gene
plotdata <- data[GLU,]
plotdata <- na.omit(plotdata)
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_20230913_Glutamine_Catabolism_Single_gene.pdf",width = 10,height = 6)
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,ncol(plotdata))]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_boxplot(aes(fill = Groupfin))+
          labs(title = coln[i]))
}
dev.off()

NEL <- c("GLS","GLS2","GLUD1","GLUD2","CPT1A","CPT1B","CPT2","CACT",
         "ACADM","ACADS","ACADVL","ACADL","ACADSB","ECH1","HADH",
         "HSD17B10","EHHADH","HSD17B4","ACAA1","ACAA2","HADHB")
plotdata <- data[NEL,]
plotdata <- na.omit(plotdata)
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_20230913_Selected_FAO&GLU_Single_gene.pdf",width = 10,height = 6)
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,ncol(plotdata))]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_boxplot(aes(fill = Groupfin))+
          labs(title = coln[i]))
}
dev.off()

PPAR <- c("PPARA","PPARD","PPARG","PPARGC1A","SREBF1","SREBF2")
plotdata <- data[PPAR,]
plotdata <- na.omit(plotdata)
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_20230913_PPAR_Single_gene.pdf",width = 10,height = 6)
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,ncol(plotdata))]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_boxplot(aes(fill = Groupfin))+
          labs(title = coln[i]))
}
dev.off()

PPARD <- c("HMGCS2","APOA1","APOA2","APOA5","PLTP","APOC3","ME1","ME3","SCD","FADS2",
           "CYP7A1","CYP8B1","NR1H3","CYP27A1","DBI","FABP1","FABP3","LPL","ACSL1",
           "ACSL3","ACSL4","ACSL5","ACSL6","ACSBG1","ACSBG2","SLC27A1","SLC27A4",
           "CD36","OLR1","EHHADH","CYP4A1","ACAA1","SCP2","ACOX1","ACOX2","ACOX3","CPT1A",
           "CPT1B","CPT1C","CPT2","ACADL","ACADM","ANGPTL4","FABP4","SORBS1","PLIN1","PLIN2",
           "PLIN4","PLIN5","ADIPOQ","MMP1","UCP1","ILK","PDPK1","UBC","PCK1","PCK2","GK",
           "GK2","AQP7","AQP7B")
plotdata <- data[PPARD,]
plotdata <- na.omit(plotdata)
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_20230913_PPAR_Downstream_Single_gene.pdf",width = 10,height = 6)
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,ncol(plotdata))]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_boxplot(aes(fill = Groupfin))+
          labs(title = coln[i]))
}
dev.off()

PPARD <- as.data.frame(t(PPARD))
PPARD <- add_column(PPARD,.before = "V1",loc = NA)
rownames(PPARD) <- "PPAR_Downstream_Genes"
write.table(PPARD,"GMTs/PPAR_Downstream_Genes.gmt",sep = "\t",row.names = T,col.names = F,quote = F)

gmtgene <- read.gmt("GMTs/PPAR_Downstream_Genes.gmt")

pdf(file = "XXY_20230913_PPAR_Downstream.pdf",height = 10,width = 14)
for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs E28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("B6I vs 28Z",Round[i],sep = " ")))
}

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("E28Z",Round[i],sep = "_"),paste("28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()


