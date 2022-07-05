library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)


load("CART_TPMvalue.RData")

IPcelllist <- read.table("GSE125881_RAW/GSM3583885_NHL-6_IP-unpaired-clonotypes.csv.gz",sep = ",",header = T)
D12celllist <- read.table("GSE125881_RAW/GSM3583888_NHL-6_d102-unpaired-clonotypes.csv.gz",sep = ",",header = T)

IPrank <- table(IPcelllist$frequency)
IPrank[length(IPrank):1] <- 1:length(IPrank)
IPcelllist$rank <- sapply(1:nrow(IPcelllist),function(x){IPrank[which(names(IPrank) == IPcelllist$frequency[x])]})

D12rank <- table(D12celllist$frequency)
D12rank[length(D12rank):1] <- 1:length(D12rank)
D12celllist$rank <- sapply(1:nrow(D12celllist),function(x){D12rank[which(names(D12rank) == D12celllist$frequency[x])]})

IPcell <- IPcelllist[which(IPcelllist$cdr3s %in% D12celllist$cdr3s),]
D12cell <- D12celllist[which(D12celllist$cdr3s %in% IPcelllist$cdr3s),]

Merge <- merge(IPcell,D12cell,by = "cdr3s" )
Merge$rank <- Merge$rank.x - Merge$rank.y
Merge <- Merge[,c("cdr3s","barcodes.x","barcodes.y","rank")]
summary <- summary(Merge$rank)
plot(Merge$rank)


Hiprolif <- Merge[which(Merge$rank > summary[5]),]
Hiprolif$barcodes <- paste(Hiprolif$barcodes.x,Hiprolif$barcodes.y,sep = ";")
Loprolif <- Merge[which(Merge$rank < summary[2]),]
Loprolif$barcodes <- paste(Loprolif$barcodes.x,Loprolif$barcodes.y,sep = ";")

HicellD12 <- str_split(Hiprolif$barcodes,";") 
HicellD12 <- unlist(HicellD12) %>% str_replace_all("-",".")
LocellD12 <- str_split(Loprolif$barcodes,";")
LocellD12 <- unlist(LocellD12) %>% str_replace_all("-",".")

Hidata <- raw[,which(colnames(raw) %in% HicellD12)]
Lodata <- raw[,which(colnames(raw) %in% LocellD12)]

