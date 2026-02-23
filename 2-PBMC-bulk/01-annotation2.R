
# set environment
rm(list = ls())
library(tidyverse)
library(data.table)
library(GEOquery)

# read exp
exp = fread(input = "./0-Rawdata/ADNI/adni.matrix.csv", data.table = F)

exp = exp[!duplicated(exp$GeneID), ]
exp = exp[!str_detect(string = exp$GeneID, pattern = " "), ]
rownames(exp) = exp$GeneID
exp = exp[,-1]

# match phenotype
ph = read.csv(file = "./0-Rawdata/ADNI/adni.ph.csv")
ph = ph[ph$DX.bl %in% c("CN", "AD"),]

rownames(ph) = ph$PTID

kk = intersect(rownames(ph), colnames(exp))

ph = ph[kk,]
exp = exp[,kk]

# save
write.csv(exp, file = "./1-Cleandata/ADNI.exp.csv")
write.csv(ph, file = "./1-Cleandata/ADNI.ph.csv")

