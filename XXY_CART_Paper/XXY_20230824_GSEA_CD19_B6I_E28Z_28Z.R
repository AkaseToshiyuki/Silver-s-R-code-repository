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

gmtgene <- bind_rows(read.gmt("GMTs/Automerge_20230102.gmt"),
                     read.gmt("GMTs/AutophagyCore_20230102.gmt"))

pdf(file = "XXY_20230824_Autophagy.pdf",height = 10,width = 14)
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

#Chol deficient
defup <- read.csv("YCS_CancerCell/YCS_BULK_20221009/up_GWRGXRO_list_1229.csv")
defup <- defup[which(defup$Gene != ""),]
defup <- as.data.frame(t(defup$Gene))
defup <- add_column(defup,.before = "V1",loc = NA)
rownames(defup) <- "GWRGXRO_Up"
write.table(defup,"GMTs/Chol_deficient_GWRGXRO_Up.gmt",sep = "\t",row.names = T,col.names = F,quote = F)

defup <- read.csv("YCS_CancerCell/YCS_BULK_20221009/down_GWRGXRO_list_1229.csv")
defup <- defup[which(defup$Gene != ""),]
defup <- as.data.frame(t(defup$Gene))
defup <- add_column(defup,.before = "V1",loc = NA)
rownames(defup) <- "GWRGXRO_Down"
write.table(defup,"GMTs/Chol_deficient_GWRGXRO_Down.gmt",sep = "\t",row.names = T,col.names = F,quote = F)


gmtgene <- bind_rows(read.gmt("GMTs/Chol_deficient_GWRGXRO_Up.gmt"),read.gmt("GMTs/Chol_deficient_GWRGXRO_Down.gmt"))

pdf(file = "XXY_20230824_Chol_Updown.pdf",height = 10,width = 14)
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


Group <- c("E28Z","B6I","28Z")
Round <- c("Ctrl","R1","R2","R3","R5")

pdf(file = "XXY_20230824_Chol_Updown_R3vsCtrl.pdf",height = 10,width = 14)
for (i in 1:3) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste(Group[i],"R3",sep = "_"),paste(Group[i],"Ctrl",sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("R3 vs Rest",Group[i],sep = " ")))
}

dev.off()

gmtgene <- read.gmt("GMTs/GESAgene20221009_YCS.gmt")
pdf(file = "XXY_20230824_Chol_Updown_Acq.pdf",height = 10,width = 14)
for (i in 1:3) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste(Group[i],"R5",sep = "_"),paste(Group[i],"Ctrl",sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:length(gseares@result$ID)),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("R5 vs Rest",Group[i],sep = " ")))
}

dev.off()