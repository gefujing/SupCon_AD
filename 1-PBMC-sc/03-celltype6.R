
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
load(file = "./3-celltype/sce.all.filt.harmony.celltype.Rdata")
load(file = "./1-qc/mycol.Rdata")
table(sce.all.filt.harmony$celltype)

# 富集分析
sce.all.filt.harmony = RunDEtest(srt = sce.all.filt.harmony, group_by = "celltype", fc.threshold = 1, only.pos = FALSE)
sce.all.filt.harmony = RunEnrichment(srt = sce.all.filt.harmony, group_by = "celltype", db = "GO_BP", species = "Homo_sapiens",
                       DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = sce.all.filt.harmony, group_by = "celltype", plot_type = "comparison", topTerm = 3)

EnrichmentPlot(
  srt = sce.all.filt.harmony, group_by = "celltype", group_use = c("Bcell", "Monocyte/Macrophage", "NKcell", "NKTcell", "ProNKcell", "Tcell"),
  plot_type = "bar"
)
ggsave(filename = "./3-celltype/12-enrichment1.pdf", width = 16.5, height = 5.5)


sce.all.filt.harmony1 = subset(sce.all.filt.harmony, Group ==  "HC")
sce.all.filt.harmony2 = subset(sce.all.filt.harmony, Group ==  "AD")
sce.all.filt.harmony1$Age = factor(sce.all.filt.harmony1$Age, levels = c("47", "59", "60", "61", "62", "64", "66", "68", "72", "73", "74", "75", "78", "79", "81", "82", "83", "84", "85", "87", "89"))
sce.all.filt.harmony2$Age = factor(sce.all.filt.harmony2$Age, levels = c("47", "59", "60", "61", "62", "64", "66", "68", "72", "73", "74", "75", "78", "79", "81", "82", "83", "84", "85", "87", "89"))
sce.all.filt.harmony2$APOE = factor(sce.all.filt.harmony2$APOE, levels = c("E3/E3", "E3/E4", "E4/E4"))


CellStatPlot(sce.all.filt.harmony2, stat.by = "celltype", group.by = "APOE", plot_type = "area", palcolor = col_vec)
ggsave(filename = "./3-celltype/13-propotion.pdf", width = 3.75, height = 3)



DEGs = sce.all.filt.harmony@tools$DEtest_celltype$AllMarkers_wilcox
DEGs = DEGs[with(DEGs, avg_log2FC > 2 & p_val_adj < 0.05), ]

# Annotate features with transcription factors and surface proteins
ht = FeatureHeatmap(
  srt = sce.all.filt.harmony, group.by = "celltype", features = DEGs$gene, feature_split = DEGs$group1,
  species = "Homo_sapiens", db = c("GO_BP", "KEGG"), anno_terms = TRUE,
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 4
)
print(ht$plot)













