library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(limma)
library(edgeR)

data <- read.table("YCS_BULK/Table ...KO mice(1).csv",sep = ",",header = T)
#表达矩阵的构建
data <- data[which(!duplicated(data$GENE)),]
rownames(data) <- data$GENE
data <- data[-21901,-1]
data <- round(data)
#分组信息录入并将分组信息加入DGE对象中
dtn <- DGEList(data)
group <- factor(rep(c("WT","KO"),each = 3),levels = c("WT","KO"))
dtn$samples$group <- group
#去除低表达基因
left <- filterByExpr(dtn)
dtn <- dtn[left, , keep.lib.sizes=FALSE]
dtn <- calcNormFactors(dtn,method = "TMM")
#计算标准化后的Log2CPM并作图
log2cpm <- cpm(dtn,log = T)
boxplot(log2cpm)
#将表达矩阵引入Limma进行分析
expr <- data
#设计实验方案
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2)))
colnames(design) <- c("WT", "KO")
rownames(design) <- colnames(expr)
#创建对比方案建立比较矩阵
contrast.matrix <- makeContrasts(WT-KO, levels=design)
#Voom处理数据
expr1 <- voom(expr,design,plot = T)
boxplot(expr1$E)
#建立线性模型
fit <- lmFit(expr1, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
#使用eBayes函数进行经验贝叶斯方法校正
fit2 <- eBayes(fit2,trend=TRUE)
#输出高变基因，BH校正Pvalue
top.table <- topTable(fit2, adjust="BH", n = Inf, coef = 1)

plotdata <- top.table
cut_off_pvalue = 0.05  #统计显著性cut off
cut_off_logFC = 1     #差异倍数值cut off
plotdata$change = ifelse(plotdata$P.Value < cut_off_pvalue & abs(plotdata$logFC) >= cut_off_logFC, 
                     ifelse(plotdata$logFC> cut_off_logFC ,'Up','Down'),
                     'Stable')
#ggplot火山图
ggplot(plotdata, aes(x = logFC, y = -log10(P.Value), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

