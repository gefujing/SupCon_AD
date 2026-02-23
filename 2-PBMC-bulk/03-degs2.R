
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(limma)
library(FactoMineR)
library(factoextra)
library(GSVA)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(ggpubr)

load(file = "./2-Batch/ad.merge.Rdata")
load(file = "./1-Cleandata/mycol.Rdata")

# geneset
cg = read.csv(file = "./3-degs/1-allcell-marker.csv")
# cg = cg[!duplicated(cg$gene),]

gs = list(
  CD16_Mono = cg$gene[cg$cluster == "CD16+Mono"][1:50],
  cDC2 = cg$gene[cg$cluster == "cDC2"][1:50],
  HLA_DPB1_CD14_Mono = cg$gene[cg$cluster == "HLA-DPB1+CD14+Mono"][1:50],
  IL1B_CD14_Mono = cg$gene[cg$cluster == "IL1B+CD14+Mono"][1:50],
  LGALS2_CD14_Mono = cg$gene[cg$cluster == "LGALS2+CD14+Mono"][1:50],
  MX1_CD14_Mono = cg$gene[cg$cluster == "MX1+CD14+Mono"][1:50],
  pDC = cg$gene[cg$cluster == "pDC"][1:50],
  S100A12_CD14_Mono = cg$gene[cg$cluster == "S100A12+CD14+Mono"][1:50]
)


# 运行gsva
X = as.matrix(exp)
gsvaPar = gsvaParam(X, gs, kcdf = "Gaussian")
gsva.es = gsva(gsvaPar)
gsva.es[1:4, 1:4]


# 可视化
identical(colnames(gsva.es), ph$id)

ha_col = HeatmapAnnotation(
  Group = ph$group,
  GSE = ph$gse,
  col = list(
    Group = c("HC" = "skyblue", "AD" = "tomato"),
    GSE = c("ADNI" = "blue", "GSE140829" = "red"),
    Score = colorRamp2(c(0, 1), c("white", "darkgreen"))
  )
)

Heatmap(gsva.es,
        top_annotation = ha_col,
        cluster_rows = T,
        cluster_columns = T,
        show_column_names = F,
        column_split = ph$GSE)


## 可视化2
pdata = as.data.frame(t(gsva.es))
pdata = pdata[,c(1,4,5,8)]
pdata$group = ph$group
pdata$gse = ph$gse

pdata = pivot_longer(pdata, cols = 1:4, names_to = "signature", values_to = "score")

p2 = ggplot(pdata, aes(x = signature, y = score, color = group))+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 1, alpha=0.5)+
  geom_boxplot(aes(fill=group), alpha=0.1)+
  scale_color_manual(values = pal_npg('nrc')(9)[c(2,1)])+
  scale_fill_manual(values = pal_npg('nrc')(9)[c(2,1)])+
  ylab("Relative score")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p2 <- p2 + stat_compare_means(aes(group = group),
                              label="p.signif",
                              show.legend = F)
p2

ggsave(p2, filename = "./3-degs/6-score.pdf", width = 4.3, height = 3.8)


