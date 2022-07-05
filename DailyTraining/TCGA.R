setwd("E:/ÎÄÏ××ÊÁÏ/°©Ö¢/TCGA data/LUAD_TGFB")
dir.create("LUAD_LDHA_all_rna")
library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.rnaseq)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(survival)
library(survminer)

#check datasets
checkTCGA("DataSets",cancerType = "LUAD",date = tail(checkTCGA("Dates"))[6])

#download datasets
downloadTCGA(cancerTypes = "LUAD",
             dataSet = "mRNAseq_Preprocess.Level_3",
             destDir = "LUAD_LDHA_all_rna",
             date = tail( checkTCGA('Dates'), 2 )[1])
Origin <- read.table("LUAD_LDHA_all_rna/gdac.broadinstitute.org_LUAD.mRNAseq_Preprocess.Level_3.2015110100.0.0/LUAD.mRNAseq_raw_counts.txt",
                     sep = "\t",header = T)

#input datasets
expr <- expressionsTCGA(LUAD.rnaseq,LUSC.rnaseq) %>% as.data.frame()
rownames(expr) <- expr[,1]
expr <- expr[,-1]
cancer <- expr[which(str_sub(rownames(expr),14,15)=="01"),]
paracancer <- expr[which(str_sub(rownames(expr),14,15)=="11"),]
index <- str_split(colnames(cancer),"\\|",simplify = T) %>% as.data.frame()
rownames(index) <- index$V2

#correct the batch effect and normalization
library(DESeq2)
group <- factor(data$dataset)
dds <- DESeqDataSetFromMatrix(apply(t(data[-20532]), 2, as.integer),DataFrame(group),~group)
dds_result <- DESeq(dds)
res <- results(dds_result,contrast = c("group","LUSC.rnaseq","LUAD.rnaseq"))
results <- res[,c(2,6)] %>% as.data.frame()
rownames(results) <- colnames(data[,-20532])
Count_new <- counts(dds_result,normalize=T)
rownames(Count_new) <- str_split(rownames(results),"\\|",simplify = T)[,2]

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
normdata$dataset <- data$dataset

#gene expression vlnplot
gene <- c("HIF1A","RAB7A","PTEN","GABARAP","KRAS","ID2","CCND1","TGFB1","SQSTM1","ID3","BMP4")
pdf("relative_fpkm_selected_genes.pdf",height = 3,width = 3)
for (i in 1:length(gene)) {
  id <- index[which(index$V1==gene[i]),"V2"]
  p <- ggplot(normdata,aes(x=`dataset`,y=log2(normdata[,id])))+
    geom_violin(aes(fill=`dataset`))+
    geom_boxplot(fill="white",width=0.2,outlier.alpha = 0)+
    theme_classic()+xlab("TCGA_dataset")+ylab("log2fpkm")+ggtitle(gene[i])+
    stat_compare_means(method = "t.test",label.x.npc = 0.25)+theme(legend.position = "none")
  print(p)
}
dev.off()

#survival curve
rownames(data) <- str_sub(rownames(data),1,12)
suv <- survivalTCGA(LUSC.clinical,extract.cols = "admin.disease_code")
rownames(suv) <- suv$bcr_patient_barcode
expression <- normdata[intersect(suv$bcr_patient_barcode,rownames(normdata)),-20532]
expression <- na.omit(t(expression))
hl <- lapply(1:nrow(expression),function(x){ifelse(expression[x,]>median(expression[x,]),"high","low")})
names(hl) <- index[rownames(expression),"V1"]
LUSC_suv <- cbind(suv[colnames(expression),],as.data.frame(hl))

pdf("LUSC_survival_curve_PTEN.pdf",height = 4.5,width = 4.5)
suvplot <- survfit(Surv(times,patient.vital_status)~ZNF90,data = LUAD_suv)
ggsurvplot(suvplot,pval = TRUE,conf.int = T,surv.median.line = "hv",
           palette = c("#FF6666","#6666FF"),ggtheme = theme_bw())
dev.off()

save(data,expr,index,fpkm,LUAD_suv,LUSC_suv,file = "R.RData")
