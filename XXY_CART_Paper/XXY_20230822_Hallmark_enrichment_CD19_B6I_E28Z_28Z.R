library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)

data <- read.csv("Others/XXY_CART/5th_z28vsB_allDEG (version 1).csv")
data <- data[,1:8]
data <- data[order(data$logFC,decreasing = T),]

B6Iup <- data$ensembl_gene_id[which(data$logFC <= (-1) & data$FDR <= 0.1)]
B6Idown <- data$ensembl_gene_id[which(data$logFC >= 1 & data$FDR <= 0.1)]

humanhallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name,ensembl_gene)

pdf(file = "XXY_20230822_B6Ivs28Z_Enrichment.pdf",width = 8,height = 12)
HM1 <- enricher(B6Iup,
                TERM2GENE = humanhallmark,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.5,
                pvalueCutoff = 0.5,
                minGSSize = 10,
                maxGSSize = 500)
dotplot(HM1,showCategory = 15) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(hjust = 1,size = 12))+
  labs(x="",y="",title = "B6Iup Hallmark")

HM2 <- enricher(B6Idown,
                TERM2GENE = humanhallmark,
                pAdjustMethod = "BH",
                qvalueCutoff = 0.5,
                pvalueCutoff = 0.5,
                minGSSize = 10,
                maxGSSize = 500)
dotplot(HM2,showCategory = 20) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(hjust = 1,size = 12))+
  labs(x="",y="",title = "B6Idown Hallmark")


KeggUp <- bitr(B6Iup,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
KeggUp <- KeggUp$ENTREZID
Kegg1 <- enrichKEGG(gene = KeggUp,organism = "hsa",pAdjustMethod = "BH",
                    pvalueCutoff = 0.5,qvalueCutoff = 0.5,
                    minGSSize = 5,maxGSSize = 500)
dotplot(Kegg1,showCategory = 15) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(hjust = 1,size = 12))+
  labs(x="",y="",title = "B6Iup Kegg")

Keggdown <- bitr(B6Idown,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
Keggdown <- Keggdown$ENTREZID
Kegg2 <- enrichKEGG(gene = Keggdown,organism = "hsa",pAdjustMethod = "BH",
                    pvalueCutoff = 0.5,qvalueCutoff = 0.5,
                    minGSSize = 5,maxGSSize = 500)
dotplot(Kegg2,showCategory = 20) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(hjust = 1,size = 12))+
  labs(x="",y="",title = "B6Idown Kegg")


GO1 <- enrichGO(gene = B6Iup,keyType = "ENSEMBL",
                OrgDb = "org.Hs.eg.db",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                minGSSize = 5,
                maxGSSize = 500)
dotplot(GO1,showCategory = 20) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(hjust = 1,size = 12))+
  labs(x="",y="",title = "B6Iup GO-BP")

GO2 <- enrichGO(gene = B6Idown,keyType = "ENSEMBL",
                OrgDb = "org.Hs.eg.db",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                minGSSize = 5,
                maxGSSize = 500)
dotplot(GO2,showCategory = 20) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(hjust = 1,size = 12))+
  labs(x="",y="",title = "B6Idown GO-BP")

dev.off()

B6Iuphm <- data$hgnc_symbol[which(data$logFC <= (-1) & data$FDR <= 0.1)]
B6Iuphm <- as.data.frame(B6Iuphm)
write.table(B6Iuphm,file = "HM.txt",quote = F,col.names = F,row.names = F)
