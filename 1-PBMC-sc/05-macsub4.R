

# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(scop)
library(patchwork)

# 设置数据
load(file = "./5-macsub/macsub.Rdata")
load(file = "./1-qc/mycol.Rdata")
table(macsub$celltype)

# 富集分析
macsub = RunDEtest(srt = macsub, group_by = "celltype", fc.threshold = 1, only.pos = FALSE)
macsub = RunEnrichment(srt = macsub, group_by = "celltype", db = "GO_BP", species = "Homo_sapiens",
                       DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = macsub, group_by = "celltype", group_use = c("LGALS2+CD14+Mono"), plot_type = "enrichmap")
ggsave(filename = "./5-macsub/12-enrichment1.pdf", width = 9, height = 6.8)

EnrichmentPlot(srt = macsub, group_by = "celltype", group_use = c("IL1B+CD14+Mono"), plot_type = "enrichmap")
ggsave(filename = "./5-macsub/12-enrichment2.pdf", width = 9, height = 6.8)

EnrichmentPlot(srt = macsub, group_by = "celltype", group_use = c("CD16+Mono"), plot_type = "enrichmap")
ggsave(filename = "./5-macsub/12-enrichment3.pdf", width = 9, height = 6.8)
ggsave(filename = "./5-macsub/12-enrichment3-2.pdf", width = 12, height = 8)

macsub1 = subset(macsub, Group ==  "HC")
macsub2 = subset(macsub, Group ==  "AD")
macsub1$Age = factor(macsub1$Age, levels = c("47", "59", "60", "61", "62", "64", "66", "68", "72", "73", "74", "75", "78", "79", "81", "82", "83", "84", "85", "87", "89"))
macsub2$Age = factor(macsub2$Age, levels = c("47", "59", "60", "61", "62", "64", "66", "68", "72", "73", "74", "75", "78", "79", "81", "82", "83", "84", "85", "87", "89"))
macsub2$APOE = factor(macsub2$APOE, levels = c("E3/E3", "E3/E4", "E4/E4"))


CellStatPlot(macsub2, stat.by = "celltype", group.by = "APOE", plot_type = "area", palcolor = col_vec)
ggsave(filename = "./5-macsub/13-propotion.pdf", width = 3.75, height = 3)


