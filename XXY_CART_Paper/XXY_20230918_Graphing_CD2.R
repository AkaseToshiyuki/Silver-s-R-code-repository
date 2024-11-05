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
library(RCurl)
library(XML)
library(ggpmisc)
library(ggpubr)
library(patchwork)
library(ggrepel)

E28z <- read.csv("Others/XXY_CART/20230404_NFAT/E_readsCount(1).csv")
N28Z <- read.csv("Others/XXY_CART/20230404_NFAT/z28_readsCount(1).csv")
B6I <- read.csv("Others/XXY_CART/20230404_NFAT/B_readsCount(1).csv")

load("GeneinfoENSG.RData")
geneinfo <- geneinfo[!duplicated(geneinfo$ensembl_gene_id),]
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

for (i in 1:5) {
  temps <- results(dds,tidy = T,contrast = c("Groupfin",paste("B6I",Round[i],sep = "_"),paste("E28Z",Round[i],sep = "_")))
  temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps <- temps[which(temps$padj <= 0.1 & abs(temps$log2FoldChange)>= 1),]
  temps <- temps[order(temps$log2FoldChange,decreasing = T),]
  assign(paste("B6I_Vs_E28Z",Round[i],sep = "_"),value = temps)
}

names <- paste("B6I_Vs_E28Z",Round,sep = "_")

#根据xpath获取节点内容：
getNodesTxt <- function(html_txt1,xpath_p){
  els1 = getNodeSet(html_txt1, xpath_p)
  #获得Node的内容，并且去除空字符：
  els1_txt <- sapply(els1,xmlValue)[!(sapply(els1,xmlValue)=="")]
  #去除\n：
  str_replace_all(els1_txt,"(\\n )+","")
}

#处理节点格式，为character且长度为0的赋值为NA：
dealNodeTxt <- function(NodeTxt){
  ifelse(is.character(NodeTxt)==T && length(NodeTxt)!=0 , NodeTxt , NA)
}


for (i in 1:5) {
  temps <- get(names[i])
  tempsup <- temps[which(temps$log2FoldChange > 0),]
  tempsdown <- temps[which(temps$log2FoldChange < 0),]
  
  genes <- as.data.frame(temps$hgnc_symbol)
  colnames(genes) <- "SYMBOL"
  genes <- bitr(genes$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  genes$NCBI_url <- paste("https://www.ncbi.nlm.nih.gov/gene/",genes$ENTREZID,sep="")
  head(genes)
  
  for(j in 1:nrow(genes)){
    # 获得网址：
    doc <- getURL(genes[j,"NCBI_url"])
    cat("成功获得网页！\t")
    # 获得网页内容
    html_txt1 = htmlParse(doc, asText = TRUE)
    # 获得Full Name:
    genes[j,"FullName"] <- dealNodeTxt(getNodesTxt(html_txt1,'//*[@id="summaryDl"]/dd[2]/text()'))
    cat("写入基因\t")
    # 获得HGNC ID:
    genes[j,"HGNC_ID"] <- str_replace_all(dealNodeTxt(getNodesTxt(html_txt1,'//*[@id="summaryDl"]/dd[3]/a')),"HGNC|:","")
    cat("写入HGNC_ID\t")
    # 获得Gene type:
    genes[j,"GeneType"] <- dealNodeTxt(getNodesTxt(html_txt1,'//*[@id="summaryDl"]/dd[5]/text()'))
    cat("写入GeneType\t")
    # 获得summary：
    genes[j,"Summary"] <- ifelse(length(getNodesTxt(html_txt1,'//*[@id="summaryDl"]/dd[10]/text()'))!=0,
                                 getNodesTxt(html_txt1,'//*[@id="summaryDl"]/dd[10]/text()'),NA)
    cat("写入Summary\n")
    print(paste("完成第",j,"个"))
    }
  genes1 <- merge(temps,genes,by.x = "hgnc_symbol",by.y = "SYMBOL")
  genes1 <- genes1[order(genes1$log2FoldChange,decreasing = T),]
  write.csv(x = genes1,file = paste(names[i],".csv",sep = ""))
}


#Enrich
B6Iup <- B6I_Vs_E28Z_R3$row[which(B6I_Vs_E28Z_R3$log2FoldChange > 0)]
B6Idown <- B6I_Vs_E28Z_R3$row[which(B6I_Vs_E28Z_R3$log2FoldChange < 0)]

humanhallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, ensembl_gene)

pdf(file = "XXY_20230917_B6IvsE28Z_Round3_Enrichment.pdf",width = 8,height = 12)
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
  labs(x="",y="",title = "B6I up Hallmark")

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
  labs(x="",y="",title = "B6I down Hallmark")


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
  labs(x="",y="",title = "B6I up Kegg")

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
  labs(x="",y="",title = "B6I down Kegg")


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
  labs(x="",y="",title = "B6I up GO-BP")

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
  labs(x="",y="",title = "B6I down GO-BP")

dev.off()


CD2up <- read.gmt("GMTs/CD2_Tcell_Upregulate.gmt")
CD2down <- read.gmt("GMTs/CD2_Tcell_downregulate.gmt")

coldata <- as.data.frame(colnames(counts))
coldata$Group <- rep(c("28Z","E28Z","B6I"),each = 15)
coldata$Round <- rep(rep(c("Ctrl","R1","R2","R3","R5"),each = 3),3)
coldata$Groupfin <- paste(coldata$Group,coldata$Round,sep = "_")

ref <- as.data.frame(rownames(counts))
ref <- bitr(ref$`rownames(counts)`,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
CD2up <- ref[which(ref$SYMBOL %in% CD2up$gene),1]
CD2down <- ref[which(ref$SYMBOL %in% CD2down$gene),1]

pdf("XXY_20230921_CD2_Score_UP&DOWN.pdf",width = 8,height = 6)

counts1 <- as.data.frame(t(counts))
counts1$uppsd <- sapply(1:nrow(counts1),function(x){sum(counts1[x,which(colnames(counts1) %in% CD2up)])/sum(counts1[x,]) * 10^4})
plotdata <- as.data.frame(counts1$uppsd)
plotdata$Groupfin <- coldata$Groupfin
plotdata <- plotdata[c(4:15,19:30,34:45),]
plotdata$Round <- str_extract(plotdata$Groupfin,pattern = "(?<=\\_).*")
plotdata$Group <- str_extract(plotdata$Groupfin,pattern = ".*(?=\\_)")
plotdata$Groupfin <- factor(x = plotdata$Groupfin,levels = c("28Z_R1","E28Z_R1","B6I_R1","28Z_R2","E28Z_R2","B6I_R2","28Z_R3","E28Z_R3","B6I_R3","28Z_R5","E28Z_R5","B6I_R5"))
plotdata$Group <- factor(x = plotdata$Group,levels = c("28Z","E28Z","B6I"))

print(ggplot(plotdata,aes(x = `Group`,y = `counts1$uppsd`))+
        geom_boxplot(aes(fill = `Group`),show.legend = T)+
        #stat_compare_means(aes(group = `Groupfin`),label = "p.format",paired = F,
        #                   method = "wilcox.test",hide.ns = F,comparisons = 
        #                     list(c("B6I_R5","E28Z_R5"),c("B6I_R5","28Z_R5")))+
        theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              #axis.text.x = element_text(angle=45,hjust = 1,size = 13),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 13),
              axis.title.y = element_text(size = 13))+
        facet_grid(. ~ Round)+
        labs(title = "", y = "CD2 Upregulated Signiture Score",x = ""))


counts1 <- as.data.frame(t(counts))
counts1$dnpsd <- sapply(1:nrow(counts1),function(x){sum(counts1[x,which(colnames(counts1) %in% CD2down)])/sum(counts1[x,]) * 10^4})
plotdata <- as.data.frame(counts1$dnpsd)
plotdata$Groupfin <- coldata$Groupfin
plotdata <- plotdata[c(4:15,19:30,34:45),]
plotdata$Round <- str_extract(plotdata$Groupfin,pattern = "(?<=\\_).*")
plotdata$Group <- str_extract(plotdata$Groupfin,pattern = ".*(?=\\_)")
plotdata$Groupfin <- factor(x = plotdata$Groupfin,levels = c("28Z_R1","E28Z_R1","B6I_R1","28Z_R2","E28Z_R2","B6I_R2","28Z_R3","E28Z_R3","B6I_R3","28Z_R5","E28Z_R5","B6I_R5"))
plotdata$Group <- factor(x = plotdata$Group,levels = c("28Z","E28Z","B6I"))

print(ggplot(plotdata,aes(x = `Group`,y = `counts1$dnpsd`))+
        geom_boxplot(aes(fill = `Group`),show.legend = T)+
        #stat_compare_means(aes(group = `Groupfin`),label = "p.format",paired = F,
        #                   method = "wilcox.test",hide.ns = F,comparisons = 
        #                     list(c("B6I_R5","E28Z_R5"),c("B6I_R5","28Z_R5")))+
        theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              #axis.text.x = element_text(angle=45,hjust = 1,size = 13),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 13),
              axis.title.y = element_text(size = 13))+
        facet_grid(. ~ Round)+
        labs(title = "", y = "CD2 Downregulated Signiture Score",x = ""))

dev.off()

#intersect
CD2up <- read.gmt("GMTs/CD2_Tcell_Upregulate.gmt")[,2]
CD2down <- read.gmt("GMTs/CD2_Tcell_downregulate.gmt")[,2]

R5BVsE_up <- intersect(B6I_Vs_E28Z_R5$hgnc_symbol[which(B6I_Vs_E28Z_R5$log2FoldChange >0)],CD2up)
R5BVsE_down <- intersect(B6I_Vs_E28Z_R5$hgnc_symbol[which(B6I_Vs_E28Z_R5$log2FoldChange <0)],CD2down)


temps <- results(dds,tidy = T,contrast = c("Groupfin","B6I_R5","28Z_R5"))
temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
temps <- temps[which(temps$padj <= 0.1 & abs(temps$log2FoldChange)>= 1),]
temps <- temps[order(temps$log2FoldChange,decreasing = T),]
tempsup <- temps[which(temps$log2FoldChange > 0),]
tempsdown <- temps[which(temps$log2FoldChange < 0),]
up_up <- as.data.frame(intersect(tempsup$hgnc_symbol,CD2up))
up_down <- as.data.frame(intersect(tempsup$hgnc_symbol,CD2down))
down_up <- as.data.frame(intersect(tempsdown$hgnc_symbol,CD2up))
down_down <- as.data.frame(intersect(tempsdown$hgnc_symbol,CD2down))

temps <- results(dds,tidy = T,contrast = c("Groupfin","B6I_R5","E28Z_R5"))
temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
temps <- temps[which(temps$padj <= 0.1 & abs(temps$log2FoldChange)>= 1),]
temps <- temps[order(temps$log2FoldChange,decreasing = T),]
tempsup <- temps[which(temps$log2FoldChange > 0),]
tempsdown <- temps[which(temps$log2FoldChange < 0),]
up_up <- as.data.frame(intersect(tempsup$hgnc_symbol,CD2up))
up_down <- as.data.frame(intersect(tempsup$hgnc_symbol,CD2down))
down_up <- as.data.frame(intersect(tempsdown$hgnc_symbol,CD2up))
down_down <- as.data.frame(intersect(tempsdown$hgnc_symbol,CD2down))

#Volcano Plot
selegenes <- c("NR4A1","SGPP2","FCER1G","RGS16","CCL22","NR4A2","EGR2","ID3")
vollist <- c("B6IE28Z","B6I28Z","E28Z28Z")
Titlelist <- c("B6I Vs E28Z 5th Round","B6I Vs 28Z 5th Round","E28Z Vs 28Z 5th Round")

B6IE28Z <- results(dds,tidy = T,contrast = c("Groupfin","B6I_R5","E28Z_R5")) %>% na.omit()
B6I28Z <- results(dds,tidy = T,contrast = c("Groupfin","B6I_R5","28Z_R5")) %>% na.omit()
E28Z28Z <- results(dds,tidy = T,contrast = c("Groupfin","E28Z_R5","28Z_R5")) %>% na.omit()

pdf("XXY_20230921_Volcano_with_label.pdf",width = 8,height = 7)
for (i in 1:length(vollist)) {
  temps <- get(vollist[i])
  temps$change <- ifelse(temps$padj < 0.1 & abs(temps$log2FoldChange) >= 1,
                           ifelse(temps$log2FoldChange > 1 ,'Up','Down'),'Stable')
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps$label <- ifelse((temps$hgnc_symbol %in% selegenes) & temps$change == "Down" , as.character(temps$hgnc_symbol),"")
  temps <- temps[which(-log10(temps$padj) <= 100 & abs(temps$log2FoldChange) <= 10),]
  
  print(ggplot(temps, aes(x = log2FoldChange,y = -log10(padj),colour=change)) +
          geom_point(alpha=0.1, size=1.5) +
          geom_point(data = temps[temps$label != "",], color = "#00468BFF",size = 2.5) +
          scale_color_manual(values=c("#00468BFF", "#d2dae2","#ED0000FF"))+
          geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
          geom_hline(yintercept = -log10(0.1),lty=4,col="black",lwd=0.8) +
          labs(x="Log2 (FoldChange)",y="-Log10 (P.adj)",title = Titlelist[i])+
          theme_bw()+
          theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank()) +
          geom_text_repel(aes(label = label),data = temps,size = 3,show.legend = F,color = "black",
                          direction="both",min.segment.length = 0.05,segment.alpha=0.6,
                          box.padding = 0.5,max.overlaps = 1400,nudge_x = 0.2,nudge_y=0.2))
}

dev.off()



#Stemness
Stem <- read.gmt("GMTs/TCell_Stem_Dysfunction_20230102.gmt")
Stem <- Stem$gene[which(Stem$term == "Naïve/Memory_like_module_100")]

temps <- results(dds,tidy = T,contrast = c("Groupfin","B6I_R5","E28Z_R5"))
temps <- temps[which(is.na(temps$log2FoldChange) == FALSE),]
temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
temps <- temps[which(temps$padj <= 0.1 & abs(temps$log2FoldChange)>= 1),]
temps <- temps[order(temps$log2FoldChange,decreasing = T),]
tempsup <- temps[which(temps$log2FoldChange > 0),]
up_up <- as.data.frame(intersect(tempsup$hgnc_symbol,Stem))

#GSEA another CarelJune
gmtgene <- read.gmt("GMTs/GESAgene20221009_YCS.gmt")
gmtgene <- gmtgene[which(gmtgene$term == "Dysfunction signature from Carl June" & gmtgene$gene != ""),]
gmtgene <- bind_rows(gmtgene,read.gmt("GMTs/Carl_June_Exhaustion.gmt"))

pdf(file = "XXY_20230921_Dysfuction.pdf",width = 12,height = 7)
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


#Volcano plot 2
selegenes <- c("NR4A1","SGPP2","FCER1G","RGS16","CCL22","NR4A2","EGR2","ID3","IL7R")
vollist <- c("B6IE28Z","B6I28Z","E28Z28Z")
Titlelist <- c("B6I Vs E28Z 5th Round","B6I Vs 28Z 5th Round","E28Z Vs 28Z 5th Round")

B6IE28Z <- results(dds,tidy = T,contrast = c("Groupfin","B6I_R5","E28Z_R5")) %>% na.omit()
B6I28Z <- results(dds,tidy = T,contrast = c("Groupfin","B6I_R5","28Z_R5")) %>% na.omit()
E28Z28Z <- results(dds,tidy = T,contrast = c("Groupfin","E28Z_R5","28Z_R5")) %>% na.omit()

pdf("XXY_20231013_Volcano_with_label_NEW.pdf",width = 8,height = 7)
for (i in 1:length(vollist)) {
  temps <- get(vollist[i])
  temps$change <- ifelse(temps$padj < 0.1 & abs(temps$log2FoldChange) >= 1,
                         ifelse(temps$log2FoldChange > 1 ,'Up','Down'),'Stable')
  temps <- merge(temps,geneinfo,by.x = "row",by.y = "ensembl_gene_id",all.x = F)
  temps$label <- ifelse((temps$hgnc_symbol %in% selegenes) & temps$change != "Stable", as.character(temps$hgnc_symbol),"")
  temps <- temps[which(-log10(temps$padj) <= 100 & abs(temps$log2FoldChange) <= 10),]
  
  print(ggplot(temps, aes(x = log2FoldChange,y = -log10(padj),colour=change)) +
          geom_point(alpha=0.1, size=1.5) +
          geom_point(data = temps[temps$label != "" & temps$change == "Down",], color = "#00468BFF",size = 2.5) +
          geom_point(data = temps[temps$label != "" & temps$change == "Up",], color = "#ED0000FF",size = 2.5) +
          scale_color_manual(values=c("#00468BFF", "#d2dae2","#ED0000FF"))+
          geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
          geom_hline(yintercept = -log10(0.1),lty=4,col="black",lwd=0.8) +
          labs(x="Log2 (FoldChange)",y="-Log10 (P.adj)",title = Titlelist[i])+
          theme_bw()+
          theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank()) +
          geom_text_repel(aes(label = label),data = temps,size = 3,show.legend = F,color = "black",
                          direction="both",min.segment.length = 0.05,segment.alpha=0.6,
                          box.padding = 0.5,max.overlaps = 1400,nudge_x = 0.2,nudge_y=0.2))
  write_csv(temps,file = paste("XXY_CART_Paper/",vollist[i],"_vol.csv",sep = ""))
}

dev.off()


##Type-2
TYP2L <- c("ENSG00000113520","ENSG00000113525","ENSG00000169194","ENSG00000107485")
counts1 <- as.data.frame(t(counts))
counts1$uppsd <- sapply(1:nrow(counts1),function(x){sum(counts1[x,which(colnames(counts1) %in% TYP2L)])/sum(counts1[x,]) * 10^4})
plotdata <- as.data.frame(counts1$uppsd)
plotdata$Groupfin <- coldata$Groupfin
plotdata <- plotdata[c(4:15,19:30,34:45),]
plotdata$Round <- str_extract(plotdata$Groupfin,pattern = "(?<=\\_).*")
plotdata$Group <- str_extract(plotdata$Groupfin,pattern = ".*(?=\\_)")
plotdata$Groupfin <- factor(x = plotdata$Groupfin,levels = c("28Z_R1","E28Z_R1","B6I_R1","28Z_R2","E28Z_R2","B6I_R2","28Z_R3","E28Z_R3","B6I_R3","28Z_R5","E28Z_R5","B6I_R5"))
plotdata$Group <- factor(x = plotdata$Group,levels = c("28Z","E28Z","B6I"))

print(ggplot(plotdata,aes(x = `Group`,y = `counts1$uppsd`))+
        geom_boxplot(aes(fill = `Group`),show.legend = T)+
        #stat_compare_means(aes(group = `Groupfin`),label = "p.format",paired = F,
        #                   method = "wilcox.test",hide.ns = F,comparisons = 
        #                     list(c("B6I_R5","E28Z_R5"),c("B6I_R5","28Z_R5")))+
        theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              #axis.text.x = element_text(angle=45,hjust = 1,size = 13),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 13),
              axis.title.y = element_text(size = 13))+
        facet_grid(. ~ Round)+
        labs(title = "", y = "Type-2 Score",x = ""))

#Type1
TYP1L <- c("ENSG00000111537","ENSG00000232810","ENSG00000164400","ENSG00000073861")
counts1 <- as.data.frame(t(counts))
counts1$uppsd <- sapply(1:nrow(counts1),function(x){sum(counts1[x,which(colnames(counts1) %in% TYP1L)])/sum(counts1[x,]) * 10^4})
plotdata <- as.data.frame(counts1$uppsd)
plotdata$Groupfin <- coldata$Groupfin
plotdata <- plotdata[c(4:15,19:30,34:45),]
plotdata$Round <- str_extract(plotdata$Groupfin,pattern = "(?<=\\_).*")
plotdata$Group <- str_extract(plotdata$Groupfin,pattern = ".*(?=\\_)")
plotdata$Groupfin <- factor(x = plotdata$Groupfin,levels = c("28Z_R1","E28Z_R1","B6I_R1","28Z_R2","E28Z_R2","B6I_R2","28Z_R3","E28Z_R3","B6I_R3","28Z_R5","E28Z_R5","B6I_R5"))
plotdata$Group <- factor(x = plotdata$Group,levels = c("28Z","E28Z","B6I"))

print(ggplot(plotdata,aes(x = `Group`,y = `counts1$uppsd`))+
        geom_boxplot(aes(fill = `Group`),show.legend = T)+
        #stat_compare_means(aes(group = `Groupfin`),label = "p.format",paired = F,
        #                   method = "wilcox.test",hide.ns = F,comparisons = 
        #                     list(c("B6I_R5","E28Z_R5"),c("B6I_R5","28Z_R5")))+
        theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              #axis.text.x = element_text(angle=45,hjust = 1,size = 13),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 13),
              axis.title.y = element_text(size = 13))+
        facet_grid(. ~ Round)+
        labs(title = "", y = "Type-1 Score",x = ""))

cytotoxic <- read.gmt("GMTs/GOBP_T_CELL_MEDIATED_CYTOTOXICITY.v2023.1.Hs.gmt")[,2]
cytotoxic <- geneinfo$ensembl_gene_id[which(geneinfo$hgnc_symbol %in% cytotoxic)]
cytotoxic <- cytotoxic[which(!duplicated(cytotoxic))]
counts1 <- as.data.frame(t(counts))
cytotoxic <- cytotoxic[which(cytotoxic %in% colnames(counts1))]
counts1$uppsd <- sapply(1:nrow(counts1),function(x){sum(counts1[x,which(colnames(counts1) %in% cytotoxic)])/sum(counts1[x,]) * 10^4})
plotdata <- as.data.frame(counts1$uppsd)
plotdata$Groupfin <- coldata$Groupfin
plotdata <- plotdata[c(4:15,19:30,34:45),]
plotdata$Round <- str_extract(plotdata$Groupfin,pattern = "(?<=\\_).*")
plotdata$Group <- str_extract(plotdata$Groupfin,pattern = ".*(?=\\_)")
plotdata$Groupfin <- factor(x = plotdata$Groupfin,levels = c("28Z_R1","E28Z_R1","B6I_R1","28Z_R2","E28Z_R2","B6I_R2","28Z_R3","E28Z_R3","B6I_R3","28Z_R5","E28Z_R5","B6I_R5"))
plotdata$Group <- factor(x = plotdata$Group,levels = c("28Z","E28Z","B6I"))

print(ggplot(plotdata,aes(x = `Group`,y = `counts1$uppsd`))+
        geom_boxplot(aes(fill = `Group`),show.legend = T)+
        #stat_compare_means(aes(group = `Groupfin`),label = "p.format",paired = F,
        #                   method = "wilcox.test",hide.ns = F,comparisons = 
        #                     list(c("B6I_R5","E28Z_R5"),c("B6I_R5","28Z_R5")))+
        theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              #axis.text.x = element_text(angle=45,hjust = 1,size = 13),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 13),
              axis.title.y = element_text(size = 13))+
        facet_grid(. ~ Round)+
        labs(title = "", y = "Cytotoxic Score",x = ""))


#Th17 
TH17 <- c("ENSG00000112115","ENSG00000112116","ENSG00000138684","ENSG00000127318","ENSG00000069667","ENSG00000143365")

counts1 <- as.data.frame(t(counts))
counts1$uppsd <- sapply(1:nrow(counts1),function(x){sum(counts1[x,which(colnames(counts1) %in% TH17)])/sum(counts1[x,]) * 10^4})
plotdata <- as.data.frame(counts1$uppsd)
plotdata$Groupfin <- coldata$Groupfin
plotdata <- plotdata[c(4:15,19:30,34:45),]
plotdata$Round <- str_extract(plotdata$Groupfin,pattern = "(?<=\\_).*")
plotdata$Group <- str_extract(plotdata$Groupfin,pattern = ".*(?=\\_)")
plotdata$Groupfin <- factor(x = plotdata$Groupfin,levels = c("28Z_R1","E28Z_R1","B6I_R1","28Z_R2","E28Z_R2","B6I_R2","28Z_R3","E28Z_R3","B6I_R3","28Z_R5","E28Z_R5","B6I_R5"))
plotdata$Group <- factor(x = plotdata$Group,levels = c("28Z","E28Z","B6I"))

print(ggplot(plotdata,aes(x = `Group`,y = `counts1$uppsd`))+
        geom_boxplot(aes(fill = `Group`),show.legend = T)+
        #stat_compare_means(aes(group = `Groupfin`),label = "p.format",paired = F,
        #                   method = "wilcox.test",hide.ns = F,comparisons = 
        #                     list(c("B6I_R5","E28Z_R5"),c("B6I_R5","28Z_R5")))+
        theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              #axis.text.x = element_text(angle=45,hjust = 1,size = 13),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 13),
              axis.title.y = element_text(size = 13))+
        facet_grid(. ~ Round)+
        labs(title = "", y = "Th-17 Score",x = ""))

#Treg
#TREG <- c("ENSG00000163235","ENSG00000105329","ENSG00000136634","ENSG00000049768","ENSG00000105246","ENSG00000168811")
TREG <- c("ENSG00000123411","ENSG00000163599","ENSG00000134460","ENSG00000049768","ENSG00000163235","ENSG00000105329","ENSG00000136634")

counts1 <- as.data.frame(t(counts))
counts1$uppsd <- sapply(1:nrow(counts1),function(x){sum(counts1[x,which(colnames(counts1) %in% TREG)])/sum(counts1[x,]) * 10^4})
plotdata <- as.data.frame(counts1$uppsd)
plotdata$Groupfin <- coldata$Groupfin
plotdata <- plotdata[c(4:15,19:30,34:45),]
plotdata$Round <- str_extract(plotdata$Groupfin,pattern = "(?<=\\_).*")
plotdata$Group <- str_extract(plotdata$Groupfin,pattern = ".*(?=\\_)")
plotdata$Groupfin <- factor(x = plotdata$Groupfin,levels = c("28Z_R1","E28Z_R1","B6I_R1","28Z_R2","E28Z_R2","B6I_R2","28Z_R3","E28Z_R3","B6I_R3","28Z_R5","E28Z_R5","B6I_R5"))
plotdata$Group <- factor(x = plotdata$Group,levels = c("28Z","E28Z","B6I"))

print(ggplot(plotdata,aes(x = `Group`,y = `counts1$uppsd`))+
        geom_boxplot(aes(fill = `Group`),show.legend = T)+
        #stat_compare_means(aes(group = `Groupfin`),label = "p.format",paired = F,
        #                   method = "wilcox.test",hide.ns = F,comparisons = 
        #                     list(c("B6I_R5","E28Z_R5"),c("B6I_R5","28Z_R5")))+
        theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              #axis.text.x = element_text(angle=45,hjust = 1,size = 13),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 13),
              axis.title.y = element_text(size = 13))+
        facet_grid(. ~ Round)+
        labs(title = "", y = "Treg Score",x = ""))

#Single Genes
##TPM
load("GeneinfoENSG.RData")
geneinfoTPMid <- geneinfo[which(duplicated(geneinfo$ensembl_gene_id)),"ensembl_gene_id"]
geneinfoTPM <- geneinfo[which(geneinfo$ensembl_gene_id %in% geneinfoTPMid),]
geneinfoTPM2 <- geneinfo[which(!(geneinfo$ensembl_gene_id %in% geneinfoTPM$ensembl_gene_id)),]
genenamedup <- geneinfoTPM[which(!duplicated(geneinfoTPM$hgnc_symbol)),"ensembl_gene_id"]

## Use MAX transcript length
pb <- txtProgressBar(min = 1,max = length(genenamedup),style = 3)
for (i in 1:length(genenamedup)) {
  geneinfoTPM[which(geneinfoTPM$ensembl_gene_id == genenamedup[i]),"transcript_length"] <- 
    max(geneinfoTPM[which(geneinfoTPM$ensembl_gene_id == genenamedup[i]),"transcript_length"])
  setTxtProgressBar(pb,i)
}
geneinfoTPM <- geneinfoTPM[which(!duplicated(geneinfoTPM$ensembl_gene_id)),]
geneinfofin <- bind_rows(geneinfoTPM,geneinfoTPM2)

# geneinfofin is the dataframe used for calculate TPM value
counts$ensembl_gene_id <- rownames(counts)
data1 <- merge(counts,geneinfofin)

pb <- txtProgressBar(min = 2,max = ncol(data1)-2,style = 3)
for (i in 2:(ncol(data1)-2)) {
  data1[,i] <- data1[,i]/data1$transcript_length
  setTxtProgressBar(pb,i)
}

colSums(data1[,2:(ncol(data1)-2)])

pb <- txtProgressBar(min = 1,max = (ncol(data1)-2),style = 3)
for (i in 2:(ncol(data1)-2)) {
  data1[,i] <- data1[,i]/sum(data1[,i])*10^6
  setTxtProgressBar(pb,i)
}

colSums(data1[,2:(ncol(data1)-2)])

#data2 <- data1[which(data1$ensembl_gene_id %in% c(TYP1L,TYP2L,TH17,TREG,"ENSG00000156127")),]
#data2 <- data1[which(data1$ensembl_gene_id %in% c("ENSG00000083457","ENSG00000138185","ENSG00000213949","ENSG00000126353")),]
data2 <- data1[which(data1$ensembl_gene_id %in% c("ENSG00000116815","ENSG00000069702",
                                                  "ENSG00000260001","ENSG00000135966",
                                                  "ENSG00000163513","ENSG00000106799",
                                                  "ENSG00000081059","ENSG00000150907")),]
rownames(data2) <- data2$hgnc_symbol
data2 <- data2[,2:46]
data2 <- as.data.frame(t(data2))


pdf("XXY_CART_Paper/XXY_20241016_Single_gene_TPM_3.pdf",width = 7,height = 6)
for (i in 1:ncol(data2)) {
  plotdata <- as.data.frame(data2[,i])
  colnames(plotdata) <- "TPM"
  plotdata$Groupfin <- coldata$Groupfin
  plotdata$Round <- str_extract(plotdata$Groupfin,pattern = "(?<=\\_).*")
  plotdata$Group <- str_extract(plotdata$Groupfin,pattern = ".*(?=\\_)")
  plotdata <- plotdata[which(plotdata$Round != "Ctrl"),]
  plotdata$Groupfin <- factor(x = plotdata$Groupfin,levels = c("28Z_R1","E28Z_R1","B6I_R1","28Z_R2","E28Z_R2","B6I_R2","28Z_R3","E28Z_R3","B6I_R3","28Z_R5","E28Z_R5","B6I_R5"))
  plotdata$Group <- factor(x = plotdata$Group,levels = c("28Z","E28Z","B6I"))
  
  print(ggplot(plotdata,aes(x = `Group`,y = `TPM`))+
          geom_boxplot(aes(fill = `Group`),show.legend = T)+
          #stat_compare_means(aes(group = `Groupfin`),label = "p.format",paired = F,
          #                   method = "wilcox.test",hide.ns = F,comparisons = 
          #                     list(c("B6I_R5","E28Z_R5"),c("B6I_R5","28Z_R5")))+
          theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                #axis.text.x = element_text(angle=45,hjust = 1,size = 13),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 13),
                axis.title.y = element_text(size = 13))+
          facet_grid(. ~ Round)+
          labs(title = colnames(data2)[i], y = "Expression Level (TPM)",x = ""))
}
dev.off()
