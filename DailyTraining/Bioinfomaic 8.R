library(e1071)
library(RTCGA)
library(RTCGA.rnaseq)
library(stringr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(survival)
library(survminer)
library(randomForest)


data <- read.csv("BRCA_clinicalMatrix(1).txt",sep = "\t",header = T)
data[data == ""] <- NA
pam <- as.data.frame(data$PAM50Call_RNAseq)
rownames(pam) <- data$sampleID
pamn <- dplyr::filter(pam,!is.na(pam$`data$PAM50Call_RNAseq`))
exp <- as.data.frame(expressionsTCGA(BRCA.rnaseq))
rownames(exp) <- exp[,1]
exp <- exp[,-1]
exp <- exp[,-20532]
same <- exp[which(str_sub(rownames(exp),1,15) %in% rownames(pamn)),]
pamn[,2] <- rownames(pamn)
pamntest <- pamn[which(rownames(pamn) %in% str_sub(rownames(same),1,15)),]
rownames(same) <- str_sub(rownames(same),1,15)
datatrain <- cbind(same[rownames(pamntest),],pamntest[,-2])
colnames(datatrain)[20532] <- "Type"
#Training data set make complete
var <- sapply(1:20531, function(x){var(datatrain[,x])})
plot(var)
data_select <- datatrain[,which(var>(5*10**9))]
train <- cbind(data_select[rownames(pamntest),],pamntest[,-2])
colnames(train) <- c("COL1A1","COL1A2","COL3A1","CPB1","FN1","SCGB2A2","Type")
attach(train)
#SVM
train$Type <- factor(train$Type)
model1 <- svm(Type ~.,data = train)
plot(model1, train, COL3A1 ~ COL1A2)
summary(model1)
#RandomForest
model2 <- randomForest(Type ~., train)
print(model2) 
plot(model2)
