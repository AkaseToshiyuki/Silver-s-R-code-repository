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

ICOSsig <- read.csv("ICOS_sig_frontiers.csv")
ICOSdown <- as.data.frame(t(as.data.frame(ICOSsig[which(ICOSsig$log2FC < -3),1])))
rownames(ICOSdown) <- "Tcell_ICOS_down"
ICOSdown <- add_column(ICOSdown,.before = "V1", loc = NA)
write.table(ICOSdown,"GMTs/ICOS_Tcell_downregulate.gmt",sep = "\t",row.names = T,col.names = F,quote = F)

ICOSup <- as.data.frame(t(as.data.frame(ICOSsig[which(ICOSsig$log2FC > 3),1])))
rownames(ICOSup) <- "Tcell_ICOS_up"
ICOSup <- add_column(ICOSup,.before = "V1", loc = NA)
write.table(ICOSup,"GMTs/ICOS_Tcell_upregulate.gmt",sep = "\t",row.names = T,col.names = F,quote = F)


CD28sig <- read.csv("CD28_sig_frontiers.csv")
CD28down <- as.data.frame(t(as.data.frame(CD28sig[which(CD28sig$log2FC < -3),1])))
rownames(CD28down) <- "Tcell_CD28_down"
CD28down <- add_column(CD28down,.before = "V1", loc = NA)
write.table(CD28down,"GMTs/CD28_Tcell_downregulate.gmt",sep = "\t",row.names = T,col.names = F,quote = F)

CD28up <- as.data.frame(t(as.data.frame(CD28sig[which(CD28sig$log2FC > 3),1])))
rownames(CD28up) <- "Tcell_CD28_up"
CD28up <- add_column(CD28up,.before = "V1", loc = NA)
write.table(CD28up,"GMTs/CD28_Tcell_upregulate.gmt",sep = "\t",row.names = T,col.names = F,quote = F)


gmtgene <- bind_rows(read.gmt("GMTs/ICOS_Tcell_downregulate.gmt"),
                     read.gmt("GMTs/ICOS_Tcell_upregulate.gmt"),
                     read.gmt("GMTs/CD28_Tcell_downregulate.gmt"),
                     read.gmt("GMTs/CD28_Tcell_upregulate.gmt"))
#gmtgene1 <- gmtgene[which(gmtgene$gene %in% test$V2),]

pdf(file = "XXY_20230908_ICOS&CD28_UpandDown.pdf")
for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:4),pvalue_table = T,color = "black",base_size = 13,
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
  print(gseaplot2(gseares,c(1:4),pvalue_table = T,color = "black",base_size = 13,
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
  print(gseaplot2(gseares,c(1:4),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()

#TCR Downstream Signal
gmtgene <- bind_rows(read.gmt("GMTs/GSE13738_RESTING_VS_TCR_ACTIVATED_CD4_TCELL_DN.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/GSE13738_RESTING_VS_TCR_ACTIVATED_CD4_TCELL_UP.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/PID_CD8_TCR_DOWNSTREAM_PATHWAY.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/PID_CD8_TCR_PATHWAY.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/BIOCARTA_CTL_PATHWAY.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/BIOCARTA_TCR_PATHWAY.v2023.1.Hs.gmt"))

pdf(file = "XXY_20230916_TCR_Downstream_UpandDown.pdf",width = 13,height = 7)
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


#Stem and function
gmtgene <- read.gmt("GMTs/TCell_Stem_Dysfunction_20230102.gmt")

pdf(file = "XXY_20230918_Stem&function.pdf",width = 12,height = 7)
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