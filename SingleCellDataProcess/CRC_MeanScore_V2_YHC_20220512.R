library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(ggpubr)


biosyngene <- c("SREBF2", "ACAT2", "HMGCS1", "HMGCR", "SQLE", "MBTPS1", "MBTPS2", "LSS", "FDPS", 
                "FDFT1", "MVD","MVK", "PMVK", "GGPS1", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", 
                "HSD17B7", "EBP", "SC5D", "DHCR7", "DHCR24","LDLR", "NPC1", "NPC2")
LXRgene <- c("ABCA1","ABCG1","NR1H3","NR1H2","MYLIP")
uptakegene <- c("LDLR", "NPC1", "NPC2")
SREBF2genes <- c("SREBF2","HMGCR","SQLE","LDLR")
proliferation <- unlist(read.table("GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION.v7.5.1.txt",sep = ","))

raw <- read.table("GSE146771_CRC.Leukocyte.Smart-seq2.TPM.txt.gz",sep = " ",header = T)
cell <- read.table("GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz",sep = "\t",header = T)
table(cell$Global_Cluster)
table(cell$Sub_Cluster)
table(cell$Tissue)
cell <- cell[which(cell$Tissue == "T" & str_sub(cell$Sub_Cluster,6,8) == "CD8"),]
table(cell$Sub_Cluster)

cluster <- split(cell$CellName,cell$Sub_Cluster)

geneused <- c(biosyngene,LXRgene)

ht12 <- t(raw[,which(colnames(raw) %in% cluster$`hT12_CD8-LEF1`)]) %>% as.data.frame()
ht13 <- t(raw[,which(colnames(raw) %in% cluster$`hT13_CD8-GPR183`)]) %>% as.data.frame()
ht14 <- t(raw[,which(colnames(raw) %in% cluster$`hT14_CD8-CX3CR1`)]) %>% as.data.frame()
ht15 <- t(raw[,which(colnames(raw) %in% cluster$`hT15_CD8-GZMK`)]) %>% as.data.frame()
ht16 <- t(raw[,which(colnames(raw) %in% cluster$`hT16_CD8-CD6`)]) %>% as.data.frame()
ht17 <- t(raw[,which(colnames(raw) %in% cluster$`hT17_CD8-CD160`)]) %>% as.data.frame()
ht18 <- t(raw[,which(colnames(raw) %in% cluster$`hT18_CD8-LAYN`)]) %>% as.data.frame()


ht12$Group <- "CD8-HT12-LEF1"
ht13$Group <- "CD8-HT13-GPR183"
ht14$Group <- "CD8-HT14-CX3CR1"
ht15$Group <- "CD8-HT15-GZMK"
ht16$Group <- "CD8-HT16-CD6"
ht17$Group <- "CD8-HT17-CD160"
ht18$Group <- "CD8-HT18-LAYN"

CD8data <- bind_rows(ht12,ht13,ht14,ht15,ht16,ht17,ht18)
metadata <- CD8data[,geneused] 

scoredata <- metadata[,geneused]
scoredata <- as.data.frame(sapply(1:ncol(scoredata),function(x){scoredata[,x]/mean(scoredata[,x])}))
colnames(scoredata) <- colnames(metadata)
rownames(scoredata) <- rownames(metadata)
scoredata$Group <- CD8data$Group
scoredata$pscore <- rowMeans(scoredata[,which(colnames(scoredata) %in% biosyngene)])
scoredata$pscore <- log2(scoredata$pscore+1)
scoredata$nscore <- rowMeans(scoredata[,which(colnames(scoredata) %in% LXRgene)])
scoredata$nscore <- log2(scoredata$nscore+1)
scoredata$Score <- scoredata$pscore - scoredata$nscore

plotdata <- scoredata[,c("Group","Score","pscore")]
plot(plotdata$Score)

ggplot(plotdata,aes(`Group`,`Score`)) + 
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic()

plotdata$IFNG <- CD8data$IFNG

ggplot(plotdata,aes(`Score`,`IFNG`)) +
  geom_point()

IFNGHi <- plotdata[which(plotdata$IFNG > 5),]
IFNGHi$Group <- "IFNGHi"
IFNGLo <- plotdata[which(plotdata$IFNG <= 5),]
IFNGLo$Group <- "IFNGLo"
IFNGbase <- rbind(IFNGHi,IFNGLo)

ggplot(IFNGbase,aes(`Group`,`Score`)) +
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic() +
  stat_compare_means(method = "t.test",label = "p.format")



progene <- colnames(CD8data[,which(colnames(CD8data) %in% proliferation)])
prolifdata <- CD8data[,progene]
prolifdata <- prolifdata[,which(colSums(prolifdata) != 0)]
prolifscore <- prolifdata

prolifscore <- as.data.frame(sapply(1:ncol(prolifscore),function(x){prolifscore[,x]/mean(prolifscore[,x])}))
colnames(prolifscore) <- colnames(prolifdata)
rownames(prolifscore) <- rownames(prolifdata)
prolifscore$Group <- CD8data$Group
prolifscore$PLIscore <- rowMeans(prolifscore[,which(colnames(prolifscore) %in% progene)])
prolifscore$PLIscore <- log2(prolifscore$PLIscore+1)
prolifscore <- prolifscore[,c("Group","PLIscore")]
plot(prolifscore$PLIscore)

plotdata2 <- cbind(prolifscore,plotdata)
plotdata2 <- plotdata2[,-which(colnames(plotdata2) == "Group")]

write.csv(plotdata2,file = "plotdata2.csv")

ggplot(plotdata2,aes(`pscore`,`PLIscore`)) +
  geom_smooth(method='lm') +
  stat_fit_glance(method = 'lm',method.args = plotdata2(formula = PLIscore ~ pscore),
                  mapping = aes(label = sprintf('R^2~"="~%.3f~~italic(P)~"="~%.2g', stat(r.squared), stat(p.value))),
                  parse = TRUE,label.x = 0.95,label.y = 0.95,family = "SH") +
  geom_point()



"
plotdata$CDK4 <- CD8data$CDK4

ggplot(plotdata,aes(`Score`,`CDK4`)) +
  geom_point()

CDK4Hi <- plotdata[which(plotdata$CDK4 > 5),]
CDK4Hi$Group <- "CDK4Hi"
CDK4Lo <- plotdata[which(plotdata$CDK4 <= 5),]
CDK4Lo$Group <- "CDK4Lo"
CDK4base <- rbind(CDK4Hi,CDK4Lo)

ggplot(CDK4base,aes(`Group`,`Score`)) +
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic() +
  stat_compare_means(method = "t.test",label = "p.format")

ggplot(CDK4base,aes(`IFNG`,`CDK4`)) +
  geom_point()
"

