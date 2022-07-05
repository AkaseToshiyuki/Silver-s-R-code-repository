setwd("/Users/silver/R Data/LIHC")
library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.rnaseq)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(survival)
library(survminer)

dir.create("CRC")
expr <- expressionsTCGA(LIHC.rnaseq) %>% as.data.frame()
rownames(expr) <- expr[,1]
expr <- expr[,-1]
cancer <- expr[which(str_sub(rownames(expr),14,14)=="0"),]
paracancer <- expr[which(str_sub(rownames(expr),14,14)!="0"),]
cancer <- cancer[which(str_sub(rownames(cancer),1,12)%in%str_sub(rownames(paracancer),1,12)),]
paracancer <- paracancer[which(str_sub(rownames(paracancer),1,12)%in%str_sub(rownames(cancer),1,12)),]
expr <- rbind(cancer,paracancer)
index <- str_split(colnames(cancer),"\\|",simplify = T) %>% as.data.frame()
rownames(index) <- index$V2
index <- index[-20532,]

data <- expr[,-20532]
expr$sample_type <- NA
expr[rownames(cancer),"sample_type"] <- "cancer"
expr[rownames(paracancer),"sample_type"] <- "paracancer"
group <- factor(expr$sample_type)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(apply(t(data), 2, as.integer),DataFrame(group),~group)
dds_result <- DESeq(dds)
res <- results(dds_result,contrast = c("group","cancer","paracancer"))
results <- res[,c(2,6)] %>% as.data.frame()
rownames(results) <- index[,2]
Count_new <- counts(dds_result,normalize=T)
rownames(Count_new) <- rownames(results)

library(biomaRt)
ensem <- useMart("ensembl",dataset = "hsapiens_gene_ensembl")
test <- getBM(attributes = c("entrezgene_id","external_gene_name","transcript_length"),mart = ensem)
test <- as.data.frame(test[,c(2,3)],test[,1])
A <- test[rownames(Count_new),]
map <- lapply(1:ncol(Count_new), function(x){A[,2]*colSums(Count_new)[x]/10^9})
names(map) <- colnames(Count_new)
map <- as.data.frame(map)
rownames(map) <- rownames(Count_new)
fpkm <- as.data.frame(Count_new/map,rownames(Count_new))
colnames(fpkm) <- str_replace_all(colnames(fpkm),"\\.","-")
normdata <- t(fpkm) %>% as.data.frame()
normdata$dataset <- group

id <- "1718"
pdf("CRC/DHCR24.pdf",height = 3,width = 3)
ggplot(normdata,aes(x=`dataset`,y=normdata[,id]))+
  geom_violin(aes(fill=`dataset`))+
  geom_boxplot(fill="white",width=0.1,outlier.alpha = 0)+
  theme_classic()+xlab("TCGA_COAD")+ylab("fpkm")+ggtitle("")+
  stat_compare_means(method = "t.test",label.x.npc = 0.25)+theme(legend.position = "none")
dev.off()
