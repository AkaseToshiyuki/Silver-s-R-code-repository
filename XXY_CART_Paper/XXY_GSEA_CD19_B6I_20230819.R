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

gmtgene <- read.gmt("GMTs/NFAT_Downstream.gmt")
CD2down <- read.csv("CD2 signature Done.csv")
CD2down <- as.data.frame(t(CD2down))
CD2down <- add_column(CD2down,.before = "V1", loc = NA)
write.table(CD2down,"GMTs/CD2_Tcell_downregulate.gmt",sep = "\t",row.names = T,col.names = F,quote = F)
gmtgene <- read.gmt("GMTs/CD2_Tcell_downregulate.gmt")

CD2Up <- read.csv("CD2 signature_Up.csv")
CD2Up <- as.data.frame(t(CD2Up))
CD2Up <- add_column(CD2Up,.before = "V1", loc = NA)
write.table(CD2Up,"GMTs/CD2_Tcell_Upregulate.gmt",sep = "\t",row.names = T,col.names = F,quote = F)
gmtgene <- read.gmt("GMTs/CD2_Tcell_Upregulate.gmt")

gmtgene <- rbind(read.gmt("GMTs/CD2_Tcell_Upregulate.gmt"),read.gmt("GMTs/CD2_Tcell_downregulate.gmt"))

pdf(file = "XXY_20230819_CD2_Merge.pdf")

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1,2),pvalue_table = T,color = c("#00468BFF","#ED0000FF"),base_size = 13,
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
  print(gseaplot2(gseares,c(1,2),pvalue_table = T,color = c("#00468BFF","#ED0000FF"),base_size = 13,
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
  print(gseaplot2(gseares,c(1,2),pvalue_table = T,color = c("#00468BFF","#ED0000FF"),base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()



pdf(file = "XXY_20230819_CD2_Merge.pdf")

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1,2),pvalue_table = T,color = "black",base_size = 13,
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
  print(gseaplot2(gseares,c(1,2),pvalue_table = T,color = "black",base_size = 13,
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
  print(gseaplot2(gseares,c(1,2),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()

#TPM
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

AP1 <- c("JUN","JUNB","JUND","FOS","FOSL1","FOSL2","BATF",
         "BATF3","ATF1","ATF2","ATF3","ATF4","ATF5","ATF6","ATF7"
         ,"BACH1","BACH2","IRF1","IRF2","IRF3","IRF4")
plotdata <- data[AP1,]
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_AP-1_parts_Expression.pdf")
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,24)]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_point()+
          labs(title = coln[i]))
}
dev.off()


#PI3K
gmtgene <- read.gmt("GMTs/HALLMARK_PI3K_AKT_MTOR_SIGNALING.v2023.1.Hs.gmt")

pdf(file = "XXY_20230820_PI3K_CD19.pdf")
for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,1,pvalue_table = T,color = "black",base_size = 13,
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
  print(gseaplot2(gseares,1,pvalue_table = T,color = "black",base_size = 13,
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
  print(gseaplot2(gseares,1,pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()

#Chol and Prolif
gmtgene <- read.gmt("GMTs/GESAgene20221009_YCS.gmt")

pdf(file = "XXY_20230820_Chol&Prolif_CD19.pdf")
for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:6),pvalue_table = T,color = "black",base_size = 13,
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
  print(gseaplot2(gseares,c(1:6),pvalue_table = T,color = "black",base_size = 13,
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
  print(gseaplot2(gseares,c(1:6),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()

#Death MTORC
gmtgene <- bind_rows(read.gmt("GMTs/MTOR_UP.N4.V1_UP.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/MTOR_UP.N4.V1_DN.v2023.1.Hs.gmt"),
                     read.gmt("GMTs/BIOCARTA_DEATH_PATHWAY.v2023.1.Hs.gmt"))

pdf(file = "XXY_20230820_Death&mTOR_CD19.pdf")
for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  gsealist <- temps$log2FoldChange
  names(gsealist) <- temps$hgnc_symbol
  gseares <- GSEA(gsealist,verbose = T,TERM2GENE = gmtgene,pvalueCutoff = 1,pAdjustMethod = "BH")
  print(gseaplot2(gseares,c(1:3),pvalue_table = T,color = "black",base_size = 13,
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
  print(gseaplot2(gseares,c(1:3),pvalue_table = T,color = "black",base_size = 13,
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
  print(gseaplot2(gseares,c(1:3),pvalue_table = T,color = "black",base_size = 13,
                  subplots = c(1,2,3),rel_heights = c(1.5,0.2,0.2),title = paste("E28Z vs 28Z",Round[i],sep = " ")))
}

dev.off()

#Chol Acq
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

acq <- read.gmt("GMTs/GESAgene20221009_YCS.gmt")
acq <- acq[which(acq$term == "Cholesterol acqisition"),2]
plotdata <- data[which(rownames(data) %in% acq),]
plotdata <- as.data.frame(t(plotdata))
plotdata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
plotdata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
plotdata$Groupfin <- paste(plotdata$Group,plotdata$Round,sep = "_")

pdf("XXY_Chol_Acq_Expression_S.pdf",width = 9)
coln <- colnames(plotdata)
for (i in 1:(ncol(plotdata)-3)) {
  plotdata1 <- plotdata[,c(i,ncol(plotdata))]
  colnames(plotdata1) <- c("Gene","Groupfin")
  print(ggplot(plotdata1,aes(x = Groupfin,y = Gene))+
          geom_point()+
          labs(title = coln[i]))
}
dev.off()

acq <- c("LDLR","SREBF2","HMGCR")


#Metabolism
load("GMTs/GMTs_Metabolism_CDpaper_YHC.RData")

pdf(file = "XXY_20230821_YHC_Metabolism_All.pdf",height = 10,width = 14)
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

#Fatty acid
gmtgene <- bind_rows(read.gmt("GMTs/WP_FATTY_ACID_BIOSYNTHESIS.v2023.1.Hs.gmt"),read.gmt("GMTs/REACTOME_FATTY_ACYL_COA_BIOSYNTHESIS.v2023.1.Hs.gmt"))

pdf(file = "XXY_20230821_FattyAcid_Biosynthesis.pdf",height = 10,width = 14)
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

#TCR
gmtgene <- bind_rows(read.gmt("GMTs/PID_CD8_TCR_DOWNSTREAM_PATHWAY.v2023.1.Hs.gmt"))

pdf(file = "XXY_20230823_TCR_DS_Signal.pdf",height = 10,width = 14)
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

#Score
TCRgene <- read.gmt("GMTs/PID_CD8_TCR_DOWNSTREAM_PATHWAY.v2023.1.Hs.gmt")
TCRgene <- merge(TCRgene,geneinfo,by.x = "gene",by.y = "hgnc_symbol",all.y = F)
TCRgene <- TCRgene$ensembl_gene_id
CD2upgene <- read.gmt("GMTs/CD2_Tcell_Upregulate.gmt")
CD2upgene <- merge(CD2upgene ,geneinfo,by.x = "gene",by.y = "hgnc_symbol",all.y = F)
CD2upgene  <- CD2upgene $ensembl_gene_id
CD2downgene <- read.gmt("GMTs/CD2_Tcell_downregulate.gmt")
CD2downgene <- merge(CD2downgene ,geneinfo,by.x = "gene",by.y = "hgnc_symbol",all.y = F)
CD2downgene  <- CD2downgene $ensembl_gene_id

plotdata <- as.data.frame(colSums(counts)/colSums(counts[which(rownames(counts)%in%CD2downgene),]))
colnames(plotdata) <- "Score"
plotdata$round <- rep(c("ctrl","R1","R2","R3","R5"),c(3,3,3,3,3))
plotdata$group <- rep(c("28Z","E28Z","B6I"),c(15,15,15))
plotdata$final <- paste(plotdata$group,plotdata$round,sep = "_")

ggplot(plotdata,aes(x = final,y = Score))+
  geom_boxplot(aes(fill = group))
