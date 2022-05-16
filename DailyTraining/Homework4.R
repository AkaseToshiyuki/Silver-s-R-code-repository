#Q1
data <- read.table("homework4-2_data.txt",sep = "\t",header = T)
patient <- data[,1]
healthy <- data[,2]
s1 <- var(patient)
s2 <- var(healthy)
n1 <- length(patient)-1
n2 <- length(healthy)-1
F <- s1/s2
F1 <- qf(0.025,n1,n2)
F2 <- qf(0.975,n1,n2)
F/F1
F/F2
pvalue <- 2*pf(F,n1,n2)
t.test(patient,healthy,var.equal = T)

#Q2
M <- c(1.8,5.8,7.1,4.6,5.5,2.4,8.3,1.2)
FM <- c(9.5,2.6,3.7,4.7,6.4,8.4,3.1,1.4)
var.test(M,FM)
t.test(M,FM,var.equal = T)
#p=0.7806,no significance
library(pwr)
sd <- (((length(M)-1)*var(M)+((length(FM)-1)*var(FM))))/((length(M)+length(FM)-2)**0.5)
dvalue <- abs(mean(M)-mean(FM))/sd

pwr.t.test(n=8,d=dvalue, sig.level = 0.05, type = "two.sample", alternative = "two.sided")
# power = 0.0510053
pwr.t.test(power = 0.80,d=dvalue, sig.level = 0.05, type = "two.sample", alternative = "two.sided")
# n = 6238.297,need 6239 samples at least

#Q3
A <- c(4.5,6.5,6,9.2,10,12,11,12,8.3)
B <- c(4,7.2,8,14,8.8,10,13,10,11.5)
shapiro.test(A)
shapiro.test(B)
var.test(A,B)
# p-value = 0.7363 var.equal=T
t.test(A,B,var.equal = T)
# p-value = 0.5764 , No significance under 95%

#Q4
pre <- c(25.46,26.10,26.14,26.68,26.96,27.04,27.15,27.64,27.11,26.37)
aft <- c(23.68,24.30,24.15,25.23,25.49,26.03,25.47,26.22,25.76,24.77)
shapiro.test(pre)
shapiro.test(aft)
var.test(pre,aft)
t.test(pre,aft,paired = T,var.equal = F)
# p-value = 2.633e-08,have significance
t.test(aft,mu=23.2,alternative = "greater",type= "one.sample")
# p-value = 2.858e-05,mean=25.11, BMI still higher than normal
dvalue <- abs(mean(aft)-23.2)/sd(aft)
pwr.t.test(n=10,d=dvalue, sig.level = 0.05, type = "one.sample", alternative = "greater")
# power = 0.4761878
beta <- 0.0587
powervalue <- 1-beta
pwr.t.test(d=dvalue,power = powervalue, sig.level = 0.001, type = "one.sample", alternative = "greater")
# n = 9.181621

#Q5
pre <- c(123,120,138,120,100,128,138,123)
aft <- c(138,116,125,136,111,132,130,110)
shapiro.test(pre)
shapiro.test(aft)
sign <- c(aft-pre)
median <- median(sign)
c <- length(which(sign > median))
u <- length(pre)/2
o <- (length(pre)/4)**0.5
z <- (c+0.5-u)/o
p_value <- 2*pnorm(z)
p_value
# p_value = 1.276326
wilcox.test(pre,aft,paired = T)
# p-value = 0.7789

#Q6
data <- read.table("homewrok4-3_data.csv",sep = ",",header = T)
rownames(data) <- data[,1]
data <- data[,-1]
result <- data.frame()
for (i in 1:nrow(data)) {
  con <- data[i,1:4]
  con <- unlist(con)
  ko <- data[i,5:8]
  ko <- unlist(ko)
  if (shapiro.test(con)$p.value > 0.05 & shapiro.test(ko)$p.value > 0.05) {
    if (var.test(con,ko)$p.value > 0.05) {result[i,1] <- t.test(con,ko,var.equal = T)$p.value
    }
    else{result[i,1] <- t.test(con,ko,var.equal = F)$p.value}
    }
  else{result[i,1] <- wilcox.test(con,ko,exact = F)$p.value}
}

result[,2] <- rownames(data)
siggene <- as.data.frame(result[which(result$V1 < 0.05),])
siggene[,2]
#significant genes
result[,3] <- p.adjust(result$V1,method = "bonferroni")
boadj <- as.data.frame(result[which(result$V3 < 0.05),])
boadj[,2]
#After bonferroni adjust significant genes
result[,4] <- p.adjust(result$V1,method = "BH")
bhadj <- as.data.frame(result[which(result$V4 < 0.05),])
bhadj[,2]
#After BH adjust significant genes

for (i in 1:nrow(data)) {
  if (condition) {
    
  }
  
}


data <- read.table("homewrok4-3_data.csv",sep = ",",header = T)
rownames(data) <- data[,1]
data <- data[,-1]
for (i in 1:nrow(data)) {
  C <- data[i,1:4] %>% unlist()
  KO <- data[i,5:8] %>% unlist()
  if(shapiro.test(C)$p.value>0.05 & shapiro.test(KO)$p.value>0.05){
    if(var.test(C,KO)$p.value>0.05){
      data$p_value[i] <- t.test(C,KO,var.equal = T)$p.value
    }else{
      data$p_value[i] <- t.test(C,KO,var.equal = F)$p.value
    }
  }else{
    data$p_value[i] <- wilcox.test(C,KO,exact = F)$p.value
  }
}
data[which(data$p_value<0.05),] %>% rownames()
data$p.adj <- p.adjust(data$p_value,method = "bonferroni")
data[which(data$p.adj<0.05),] %>% rownames()
data$p.adj <- p.adjust(data$p_value,method = "BH")
data[which(data$p.adj<0.05),] %>% rownames()

