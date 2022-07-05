library(dplyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(stringr)


biosyngene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
uptakegene <- c("LDLR", "NPC1", "NPC2")
SREBF2genes <- c("SREBF2","HMGCR","SQLE")

patient <- read.table("GSE120575_patient_ID_single_cells.txt",header = F)
header <- read.table("GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt",sep = "\t",header = FALSE,nrows = 2)
test1 <- read.table("GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt",sep = "\t",header = FALSE,skip = 2) #read the dataset

colnames(test1) <- header[1,]
colnames(test1)[1] <- "gene"
rownames(test1) <- test1$gene
test1 <- test1[,-1]
test1[1:5,1:5]

responder <- patient[which(patient$V10=="Responder" & str_detect(patient$V9,"Pre")),]
nonresponder <- patient[which(patient$V10=="Non-responder" & str_detect(patient$V9,"Pre")),]
relist <- as.data.frame(aggregate(responder$V9,list(responder$V9),length)) #or use table($) to get the list of post_
relist <- unlist(relist[,-2])
nonrelist <- as.data.frame(aggregate(nonresponder$V9,list(nonresponder$V9),length)) 
nonrelist <- unlist(nonrelist[,-2])

header <- header[,-1]
recandidate <- header[,header[2,] %in% relist] 
recellid<- unlist(recandidate[1,])
nrecandidate <- header[,header[2,] %in% nonrelist] 
nrecellid<- unlist(nrecandidate[1,])

resdata <- test1[,colnames(test1) %in% recellid]
nresdata <- test1[,colnames(test1) %in% nrecellid]

usedgene1 <- SREBF2genes
usedgene2 <- LXRgene

resdt1 <- resdata[usedgene1,]
nresdt1 <- nresdata[usedgene1,]
resdt2 <- resdata[usedgene2,]
nresdt2 <- nresdata[usedgene2,]
resdt2 <- -(resdt2)
nresdt2 <- -(nresdt2)

resdt <- rbind(resdt1,resdt2)
nresdt <- rbind(nresdt1,nresdt2)

resdt <- t(resdt) %>% as.data.frame()
nresdt <- t(nresdt) %>% as.data.frame()
resdt$Score <- rowMeans(resdt)
nresdt$Score <- rowMeans(nresdt)
resdt$Group <- "Responder"
nresdt$Group <- "Non-Responder"


metadata <- bind_rows(resdt,nresdt)
plotdata <- metadata[,c("Score","Group")]


ggplot(plotdata,aes(`Group`,`Score`)) + 
  labs(x = "Group", y = "Score", title = "Score") +
  geom_violin(aes(fill = `Score`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20")+
  theme_classic()
