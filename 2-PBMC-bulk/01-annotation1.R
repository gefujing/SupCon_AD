
# set environment
rm(list = ls())
library(tidyverse)
library(data.table)
library(GEOquery)

# read exp
exp = fread(input = "./0-Rawdata/GSE140829.matrix.csv", data.table = F)

rownames(exp) = exp$ID_REF
exp = exp[,-1]

# get annotation
gpl = getGEO("GPL15988", AnnotGPL = TRUE)
annotation = Table(gpl)

annotation = annotation[!duplicated(annotation$`GENE SYMBOL`),]
annotation = na.omit(annotation)
annotation = annotation[!str_detect(annotation$`GENE SYMBOL`, pattern = "-"),]

# match
rownames(annotation) = annotation$ID
kk = intersect(rownames(exp), rownames(annotation))

exp = exp[kk,]
annotation = annotation[kk,]
rownames(exp) = annotation$`GENE SYMBOL`

# match phenotype
ph = read.csv(file = "./0-Rawdata/GSE140829.ph.csv")
ph = ph[ph$diagnosis %in% c("Control", "AD"),]

rownames(ph) = ph$geo_accession

kk = intersect(rownames(ph), colnames(exp))

ph = ph[kk,]
exp = exp[,kk]

# save
write.csv(exp, file = "./1-Cleandata/GSE140829.exp.csv")
write.csv(ph, file = "./1-Cleandata/GSE140829.ph.csv")

