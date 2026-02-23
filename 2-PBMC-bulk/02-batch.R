
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(DESeq2)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(AnnoProbe)
library(tinyarray)
library(data.table)
library(tinyarray)
library(limma)
load(file = "./1-Cleandata/mycol.Rdata")

# 数据集
exp.GSE140829 = fread(input = "./1-Cleandata/GSE140829.exp.csv", data.table = F)
ph.GSE140829 = read.csv(file = "./1-Cleandata/GSE140829.ph.csv")

exp.ADNI = fread(input = "./1-Cleandata/ADNI.exp.csv", data.table = F)
ph.ADNI = read.csv(file = "./1-Cleandata/ADNI.ph.csv")

# 取交集
rownames(exp.GSE140829) = exp.GSE140829$V1
exp.GSE140829 = exp.GSE140829[,-1]

rownames(exp.ADNI) = exp.ADNI$V1
exp.ADNI = exp.ADNI[,-1]

ph.GSE140829 = ph.GSE140829[,c(1,4)]
ph.GSE140829$gse = "GSE140829"
ph.ADNI = ph.ADNI[,c(1,10)]
ph.ADNI$gse = "ADNI"

# 表型数据合并
colnames(ph.GSE140829) = c("id", "group", "gse")
colnames(ph.ADNI) = c("id", "group", "gse")

ph = rbind(ph.GSE140829, ph.ADNI)
table(ph$group)
table(ph$gse)

ph$group = ifelse(ph$group %in% c("CN", "Control"), "HC", "AD")
ph$group = factor(ph$group, levels = c("HC", "AD"))

# 表达矩阵合并
mrna = intersect(rownames(exp.GSE140829), rownames(exp.ADNI))

exp.GSE140829 = exp.GSE140829[mrna,]
exp.ADNI = exp.ADNI[mrna,]

exp = cbind(exp.GSE140829, exp.ADNI)

# 删除冗余数据
rm(list = c("exp.GSE140829", "exp.ADNI", "mrna",         
            "ph.GSE140829", "ph.ADNI"))


# 分布图
data = pivot_longer(data = exp, cols = 1:757)
data$datasets = ifelse(data$name %in% ph$id[ph$gse == "GSE140829"], "GSE140829", "ADNI")

ggplot(data = data, mapping = aes(x = name, y = value, color = datasets))+
  geom_violin()+
  theme_bw()
ggsave(filename = "./2-Batch/1-human-ad-before-boxplot.pdf", width = 30, height = 5)


# PCA图
grouplist = ph$gse
grouplist = factor(grouplist, levels = unique(grouplist))
batch = ph$gse

# PCA
dat = as.data.frame(t(exp))
dat.pca = PCA(dat, graph = FALSE)
pca_plot = fviz_pca_ind(dat.pca,
                        geom.ind = "point", # show points only (nbut not "text")
                        col.ind = batch, # color by groups
                        palette = col_vec[1:2],
                        addEllipses = TRUE, # Concentration ellipses
                        legend.title = "Datasets")

ggsave(pca_plot, filename = "./2-Batch/1-human-ad-before.pdf", width = 4.5, height = 3)

# 去批次效应
grouplist = ph$group
grouplist = factor(grouplist, levels = c("HC", "AD"))

design = model.matrix(~0+grouplist)
colnames(design) = levels(grouplist)
rownames(design) = colnames(exp)

batch = ph$gse
exp2 = removeBatchEffect(x = exp, batch = batch, design = design)

# PCA2
dat = as.data.frame(t(exp2))
dat.pca = PCA(dat, graph = FALSE)
pca_plot = fviz_pca_ind(dat.pca,
                        geom.ind = "point", # show points only (nbut not "text")
                        col.ind = batch, # color by groups
                        palette = col_vec[1:2],
                        addEllipses = TRUE, # Concentration ellipses
                        legend.title = "Datasets")

ggsave(pca_plot, filename = "./2-Batch/2-human-ad-after.pdf", width = 4.5, height = 3)

exp = as.data.frame(exp2)
ph = ph

# 保存数据
save(exp, ph, file = "./2-Batch/ad.merge.Rdata")

