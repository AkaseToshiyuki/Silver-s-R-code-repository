library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggpmisc)

CRlist <- c("ac01","ac05","ac07","ac08","ac09","ac10","ac12","ac14","ac16")
PDlist <- c("ac02","ac03","ac04","ac11","ac13","ac15","ac17","ac18","ac19","ac20","ac21","ac22","ac23","ac24")
plotdata <- read.table("CART_NM_plotdata.csv",sep = ",",header = T)


for (i in 1:nrow(plotdata)) {
  if (plotdata$Sample[i] %in% CRlist) {
    plotdata$Group[i] <- "CR"}
  else{
    plotdata$Group[i] <- "PD"}  
}

write.csv(plotdata,file = "plotdt_fin_20220517_CD8_scale.csv")


ggplot(plotdata,aes(`Group`,`pscoreM`)) +
  geom_violin(aes(fill = `Group`)) +
  geom_boxplot(width = 0.1,fill = "white",outlier.alpha = 0) +
  #geom_jitter(size = 0.1, width = 0.1,shape = 16,alpha = 0.5,color = "grey20") +
  theme_classic() +
  stat_compare_means(method = "t.test",label = "p.format")
