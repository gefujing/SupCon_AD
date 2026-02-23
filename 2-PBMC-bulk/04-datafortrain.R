

# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(limma)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(AnnoProbe)
library(tinyarray)
library(data.table)
library(ggrepel)
library(eulerr)
library(ggvenn)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(aPEAR)
library(GOplot)
library(GseaVis)
library(eulerr)

load(file = "./2-Batch/ad.merge.Rdata")
load(file = "./1-Cleandata/mycol.Rdata")

# 导入差异基因
deg1 = read.csv(file = "./4-datafortrain/14-DEGs/3-DEGs_S100A12_CD14_Mono_AD_vs_HC.csv")
deg1 = deg1$gene[deg1$direction != "ns"]

deg2 = read.csv(file = "./4-datafortrain/14-DEGs/3-DEGs_MX1_CD14_Mono_AD_vs_HC.csv")
deg2 = deg2$gene[deg2$direction != "ns"]

deg3 = read.csv(file = "./4-datafortrain/14-DEGs/3-DEGs_LGALS2_CD14_Mono_AD_vs_HC.csv")
deg3 = deg3$gene[deg3$direction != "ns"]

deg4 = read.csv(file = "./4-datafortrain/14-DEGs/3-DEGs_IL1B_CD14_Mono_AD_vs_HC.csv")
deg4 = deg4$gene[deg4$direction != "ns"]

deg5 = read.csv(file = "./4-datafortrain/14-DEGs/3-DEGs_HLA_DPB1_CD14_Mono_AD_vs_HC.csv")
deg5 = deg5$gene[deg5$direction != "ns"]

deg6 = read.csv(file = "./4-datafortrain/14-DEGs/3-DEGs_CD16_Mono_AD_vs_HC.csv")
deg6 = deg6$gene[deg6$direction != "ns"]

scdeg = unique(c(deg1, deg2, deg3, deg4, deg5, deg6))

scdeg = read.csv(file = "./4-datafortrain/1-allcell-marker.csv")
scdeg = scdeg$X[scdeg$cluster == "Monocyte/Macrophage"]


addeg = read.csv(file = "./3-degs/1-advshc.degs.csv")
addeg = addeg$X[addeg$change != "NOT"]

cg = unique(c(addeg, scdeg))
cg = intersect(addeg, scdeg)

exp= exp[cg,]

write.csv(exp, file = "./4-datafortrain/exp_mono.csv")
write.csv(ph, file = "./4-datafortrain//ph_mono.csv")


vd = venneuler(c(""=0.3, ""=0.3, ""=0.2))
plot(vd)

vd <- euler(c(AD_DEGs = 2201, Myeloid_Markers = 4322, "AD_DEGs&Myeloid_Markers" = 421))
plot(vd,
     fills = list(fill = c("#fbb4ae", "#b3cde3"), alpha = 0.6),
     labels = list(col = "black", font = 4), 
     edges = FALSE,
     quantities = TRUE)








