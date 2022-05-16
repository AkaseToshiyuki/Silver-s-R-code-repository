library(RTCGA)
library(RTCGA.rnaseq)
library(factoextra)
library(ggplot2)
library(stringr)

exp <- expressionsTCGA(BRCA.rnaseq) %>% as.data.frame()
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp <- exp[,-20532]
exp <- scale(exp)

expr <- exp[,which(colSums(exp)!=0)]

fviz_nbclust(expr,FUNcluster = kmeans,method = "wss")

set.seed(123)
km.res <- kmeans(expr,centers = 15,nstart = 25)
print(km.res)
fviz_cluster(km.res,expr,repel = T)

re <- lapply(1:ncol(exp), function(x){
  tmp <- exp[,x]
  res <- sapply(1:nrow(exp),function(y){is.numeric(tmp[y])})
  return(res)
  })
head(re)
is.numeric(exp[,1])
