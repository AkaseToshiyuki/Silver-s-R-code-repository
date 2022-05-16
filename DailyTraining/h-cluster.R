library(RTCGA)
library(RTCGA.rnaseq)
library(ggplot2)
library(stringr)

exp <- expressionsTCGA(BRCA.rnaseq) %>% as.data.frame()
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp <- exp[,-20532]
expr <- exp[,which(colSums(exp)!=0)]
hclust <- hclust(dist(expr))
plot(hclust,cex = 0.001)
dev.off()
