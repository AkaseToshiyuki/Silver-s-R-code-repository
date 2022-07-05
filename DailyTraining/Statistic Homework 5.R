library(dplyr)
library(stringr)
library(data.table)

#Q1(1)
data <- read.table("/RData/Data_H5/yield.csv",header = T,sep = ",", 
                   colClasses = c("factor","numeric") )
data <- transform(data, fertilizer= factor(fertilizer,levels = c("1","2","3")))
attach(data)
str(data)
grandmean <- mean(yield)
grandmean
Smeans <- aggregate(yield, by = list(fertilizer),FUN = mean)
str(Smeans)
Svars <- aggregate(yield, by = list(fertilizer),FUN = var)
Svars
Slens <- aggregate(yield, by = list(fertilizer),FUN = length)
Slens
within_SS <- sum((Slens$x-1)*Svars$x)
within_SS
total_SS <- sum((yield-grandmean)^2)
total_SS
between_SS <- total_SS - within_SS
between_SS
df_between <- length(levels(fertilizer))-1
df_between
df_within <- length(yield) - length(levels(fertilizer))
df_within
MS_between <- between_SS/df_between
MS_between
MS_within <- within_SS/df_within
MS_within
F_ratio <- MS_between/MS_within
F_ratio
p_value <- 1-pf(F_ratio,df_between,df_within)
p_value
#Q1(2)
aov.re <- aov(data=data,yield~fertilizer)
summary(aov.re)
TukeyHSD(aov.re)
"
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = yield ~ fertilizer, data = data)

$fertilizer
         diff         lwr       upr     p adj
2-1 0.1761687 -0.19371896 0.5460564 0.4954705
3-1 0.5991256  0.22923789 0.9690133 0.0006125
3-2 0.4229569  0.05306916 0.7928445 0.0208735
"
#Q1(3)
boxplot(yield~fertilizer,data= data)

#Q2(1)
data <- read.table("/RData/Data_H5/data5_4.csv",header = T,sep = ",")
dataA <- as.data.frame(data[which(data$A == "A1"),3])
colnames(dataA) <- "A1"
dataA$A2 <- data[which(data$A == "A2"),3]
dataA$A3 <- data[which(data$A == "A3"),3]
shapiro.test(dataA$A1)
shapiro.test(dataA$A2)
shapiro.test(dataA$A3)
bartlett.test(Weight ~ A, data = data)
#p_value= 0.5676 0.1989 0.6281 all greater than 0.05 
#bartlett test p_value 0.9053 so variance same
dataB <- as.data.frame(data[which(data$B == "B1"),3])
colnames(dataB) <- "B1"
dataB$B2 <- data[which(data$B == "B2"),3]
shapiro.test(dataB$B1)
shapiro.test(dataB$B2)
bartlett.test(Weight~B, data = data)
#p_value= 0.1718 0.4853 all greater than 0.05 
#bartlett test p_value 0.6471 so variance same
#So all accorring to A and B vector is Normal Distribution,could fit variance test 
#Q2(2)
aov.re <- aov(Weight ~ A*B, data =  data)
summary(aov.re)
TukeyHSD(aov.re)
"
            Df Sum Sq Mean Sq F value   Pr(>F)    
A            2   9080    4540  56.809 5.22e-14 ***
B            1    127     127   1.589    0.213    
A:B          2     17       9   0.108    0.897    
Residuals   54   4316      80                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Weight ~ A * B, data = data)

$A
       diff       lwr       upr     p adj
A2-A1 25.39 18.576905 32.203095 0.0000000
A3-A1 26.75 19.936905 33.563095 0.0000000
A3-A2  1.36 -5.453095  8.173095 0.8805321

$B
      diff       lwr      upr     p adj
B2-B1 2.91 -1.717782 7.537782 0.2128406

$`A:B`
              diff        lwr       upr     p adj
A2:B1-A1:B1  25.09  13.277922  36.90208 0.0000009
A3:B1-A1:B1  25.49  13.677922  37.30208 0.0000006
A1:B2-A1:B1   1.87  -9.942078  13.68208 0.9970585
A2:B2-A1:B1  27.56  15.747922  39.37208 0.0000001
A3:B2-A1:B1  29.88  18.067922  41.69208 0.0000000
A3:B1-A2:B1   0.40 -11.412078  12.21208 0.9999985
A1:B2-A2:B1 -23.22 -35.032078 -11.40792 0.0000050
A2:B2-A2:B1   2.47  -9.342078  14.28208 0.9892562
A3:B2-A2:B1   4.79  -7.022078  16.60208 0.8358408
A1:B2-A3:B1 -23.62 -35.432078 -11.80792 0.0000035
A2:B2-A3:B1   2.07  -9.742078  13.88208 0.9952518
A3:B2-A3:B1   4.39  -7.422078  16.20208 0.8800277
A2:B2-A1:B2  25.69  13.877922  37.50208 0.0000005
A3:B2-A1:B2  28.01  16.197922  39.82208 0.0000001
A3:B2-A2:B2   2.32  -9.492078  14.13208 0.9919360
"

#Q3(1)
data <- read.table("Data_H5/time.tsv",header = T)
boxplot(data)
summary(data)
#Q3(3)
shapiro.test(data$control)
shapiro.test(data$low)
shapiro.test(data$middle)
shapiro.test(data$high)
data <- read.table("Data_H5/time.tsv",header = T)
data <- melt(data)
bartlett.test(value~variable, data = data)
# shapiro.test p_value = 0.7874 0.3671 0.4841 0.6357 all greater than 0.05
# bartlett test p-value = 0.6329, no variance according to medicine
#Q3(4)
aov.re <- aov(value ~ variable,data = data)
summary(aov.re)
"
            Df Sum Sq Mean Sq F value   Pr(>F)    
variable     3 1205.4   401.8   26.09 1.08e-10 ***
Residuals   56  862.4    15.4 
"
TukeyHSD(aov.re)
"
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = value ~ variable, data = data)

$variable
                    diff        lwr       upr     p adj
low-control     1.921333 -1.8729069  5.715574 0.5412872
middle-control  6.824667  3.0304265 10.618907 0.0000805
high-control   11.524667  7.7304265 15.318907 0.0000000
middle-low      4.903333  1.1090931  8.697574 0.0062459
high-low        9.603333  5.8090931 13.397574 0.0000001
high-middle     4.700000  0.9057598  8.494240 0.0094149
"

#Q4
data <- read.table("Data_H5/lung_cancer.tsv",header = T)
aggregate(data$exp_A,list(data$hospital), function(x){shapiro.test(x)$p.value})
aggregate(data$exp_A,list(data$smoke), function(x){shapiro.test(x)$p.value})
aggregate(data$exp_A,list(data$drink), function(x){shapiro.test(x)$p.value})
aggregate(data$exp_A,list(data$gender), function(x){shapiro.test(x)$p.value})
aggregate(data$exp_A,list(data$status), function(x){shapiro.test(x)$p.value})
aov.re <- aov(data=data,formula = exp_A~hospital*gender*smoke*drink*status)
summary(aov.re)
aov.fin <- aov(data=data,formula = exp_A~hospital*smoke)
summary(aov.fin)
"
                Df Sum Sq Mean Sq F value   Pr(>F)    
hospital         1  29236   29236  40.180 1.56e-09 ***
smoke            1   9135    9135  12.554 0.000494 ***
hospital:smoke   1   2691    2691   3.698 0.055923 .  
Residuals      196 142617     728                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
"


data1 <- str_replace_all(data$hospital,c("A"="1","B"="2","C"="3"))
#如果分组值不是数字，替换字符为数字
