library(dplyr)
library(stringr)

biosyngene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24","LDLR", "NPC1", "NPC2")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
uptakegene <- c("LDLR", "NPC1", "NPC2")
SREBF2genes <- c("SREBF2","HMGCR","SQLE","LDLR")
functiongene <- c("IFNG", "IL10", "GZMA", "GZMB", "PRF1","TNFSF10", "FASLG", "TNF", "IL2","CD44", "IL2RA")
proliferation <- unlist(read.table("GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION.v7.5.1.txt",sep = ","))

#upper just run once

load("ac01_TPM.RData")
raw[1:10,1:10]

geneused <- c(biosyngene,LXRgene)

scaledata <- as.data.frame(t(raw))
scaledata[1:10,1:10]

plot(scaledata$CD8A)
rule <- summary(scaledata$CD8A)
scaledata <- scaledata[which(scaledata$CD8A > rule[4]),]

scaledata[1:5,biosyngene]
scaledata[1:5,LXRgene]

scaledata$pscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% biosyngene)])
scaledata$nscore <- rowMeans(scaledata[,which(colnames(scaledata) %in% LXRgene)])
scaledata$score <- scaledata$pscore - scaledata$nscore

plotdata <- scaledata[,c("score","pscore",'nscore')]
plot(plotdata$pscoreM)
summary(plotdata)
write.csv(plotdata,file = "ac01_plotdata.csv")
