p0 <- 230/500
p1 <- 0.5
zeta <- (p0-p1)/((p1*(1-p1)/500)**0.5)
qnorm(0.025,0,1)

sapply(1:5,function(x){sd(unlist(rawdata[x,]))})


samp <- c("159","280","101","212","224",
          "379","179","264","222","362",
          "168","250","149","260","485","170")
p <- length(samp[which(samp>255)])/length(samp)
q <- length(samp[which(samp<=255)])/length(samp)

qbinom(0.95,16,p)>=length(samp[which(samp>255)])

speed<-c(914, 920, 910, 934, 953, 940, 912, 924, 930)
t.test(speed, mu=950, conf.level = 0.99)

# p值=0.001<0.01, 不接收原假设
# 初速度有明显变化，且又有926.3 <950
# 可以认为弹药初速度有显著降低

binom.test(230,500,p=0.5)


x<-c(0.140,0.138,0.143,0.142,0.144,0.137)
y<-c(0.135,0.140,0.142,0.136,0.138,0.130)
var.test(x, y, conf.level = 0.99)

t.test(x,y,var.equal = TRUE)
