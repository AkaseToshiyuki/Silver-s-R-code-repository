library(Rtsne)
library(ggplot2)
#数据读取
Data <- read.csv("Data.csv")
#对数据进行初步分析
str(Data)
dim(Data)
head(Data)
#对数据进行抽样
data_sample <- Data[sample(1:nrow(Data), size = 2000),]
#进行模型训练
tsne <- Rtsne(as.matrix(Data), check_duplicates = FALSE, pca = TRUE,pca_center=TRUE,
              pca_scale=TRUE,perplexity=30, theta=0.3, dims=2,verbose=TRUE,max_iter=2000)
#将降维后需要的结果转为数据框形式
#绘制图形
aftersne = as.data.frame(tsne$Y)
colnames(aftersne) = c("tSNE1","tSNE2")
ggplot(aftersne,aes(tSNE1,tSNE2),color=rainbow()) +
  geom_point()
