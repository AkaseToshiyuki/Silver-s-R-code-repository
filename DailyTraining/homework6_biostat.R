library(stringr)
library(dplyr)
library(ISLR)
library(caTools)
library(car)
library(ggplot2)

#Q1
data <- read.table("homework6_data/homework6_1.csv",sep = ",",header = T)
func <- lm(data = data,formula = Sales~Population+Income)
summary(func)
predict(func,newdata = data.frame(Population = 200,Income = 3000),interval = "prediction",level = 0.95)

#Q2
data <- ISLR::Default
summary(data)

set.seed(123)
data$split <- sample.split(data,SplitRatio = 1/3 )
train_data <- subset(data,split == F)
sample_data <- subset(data,split == T)
loglm <- glm(data = train_data,formula = default~student+balance+income,family = binomial(link = "logit"))
summary(loglm)

loglmr <- glm(data = train_data,formula = default~balance+student,family = binomial(link = "logit"))
summary(loglmr)
aov <- anova(loglm,loglmr,test = "Chisq")
aov

sample_data$prob <- predict(loglmr,newdata = sample_data, type = "response")
sample_data$predict <- NA
sample_data[which(sample_data$prob > 0.5),"predict"] <- "Yes"
sample_data[which(sample_data$prob < 0.5),"predict"] <- "No"
sample_data$same <- NA
sample_data[which(sample_data$predict == sample_data$default),"same"] <- "Yes"
sample_data[which(sample_data$predict != sample_data$default),"same"] <- "No"
accuracy <- length(which(sample_data$same == "Yes"))/nrow(sample_data)
accuracy

#Q3
data <- read.table("homework6_data/Salary_Data.csv",sep = ",",header = T)
data$split <- sample.split(data,SplitRatio = 1/3 )
train_data <- subset(data,split == F)
test_data <- subset(data,split == T)

func <- lm(data = train_data,formula = Salary~YearsExperience)
summary(func)

pred <- as.data.frame(predict(func,newdata = test_data,interval = "prediction",level = 0.95))
pred

scatterplot(Salary~YearsExperience,data = test_data,regLine = F,smooth = F)
abline(func,col = "green")

RMSE <- test_data
RMSE$predict <- pred$fit
sqrt(mean((RMSE$Salary - RMSE$predict)^2))

#Q4
data <- cars
data[1:6,]
scatterplot(dist~speed,data = data)

ggplot(data,aes(`speed`)) + geom_density(color = "blue",fill = "blue",alpha = 0.4)
shapiro.test(data$speed)
ggplot(data,aes(`dist`)) + geom_density(color = "red",fill = "red",alpha = 0.4)
shapiro.test(data$dist)

func <- lm(data = data,formula = dist~speed)
summary(func)
plot(func)

#Q5
data <- read.table("homework6_data/lung_cancer.txt",sep = "\t",header = T)
data$exp_A <- data$exp_A/data$exp_actin
func <- glm(data = data,formula = status~exp_A+gender+smoke+drink+age,family = binomial(link = "logit"))
summary(func)
funcnew <- glm(data = data,formula = status~exp_A+smoke,family = binomial(link = "logit"))
summary(funcnew)
