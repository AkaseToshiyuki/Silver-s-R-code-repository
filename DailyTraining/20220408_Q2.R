getwd()
setwd("D:/RData") #set work space

rawdata <- read.table("homework3.1_data.txt") #read table as dataframe and name as rawdata
rownames(rawdata) 
colnames(rawdata)

nrow(rawdata) #number of row
ncol(rawdata) #number of column
head(rawdata,3) #show the top 5 row (or use rawdata[1:3,])
top5row <- head(rawdata,5) #extract top 5 row

rmean <- vector() #creat vector
rmed <- vector()
rvar <- vector()
for (i in 1:nrow(top5row)) #repeat until get all rows
{
  rmean = c(rmean,mean(as.double(top5row[i,])))
  rmed = c(rmed,median(as.double(top5row[i,])))
  rvar = c(rvar,var(as.double(top5row[i,])))
}
top5result <- cbind(rmean,rmed,rvar) #combine the data together
rownames(top5result) <- rownames(top5row)
colnames(top5result) <- c("平均值","中位数", "方差")
top5result

A0SH <-as.data.frame(rawdata[,"A0SH"])
rownames(A0SH) <- rownames(rawdata)
noexpressgenes <- as.data.frame(rownames(rawdata[which(rawdata$A0SH==0),]))
noexpressgenes

sample10<- rawdata[sample(nrow(rawdata),10),] #Take 10 random samples
sample10 <- t(sample10) #exchange row and column
boxplot(sample10) #print boxplot graph

