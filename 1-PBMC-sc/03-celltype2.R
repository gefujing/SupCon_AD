
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(ggheatmap)
library(ggsci)
library(Nebulosa)
library(RColorBrewer)
library(patchwork)
library(scRNAtoolVis)

# 设置数据
load(file = "./3-celltype/sce.all.filt.harmony.celltype.Rdata")
load(file = "./1-qc/mycol.Rdata")

# 标记基因
Idents(sce.all.filt.harmony) = sce.all.filt.harmony@meta.data$celltype
sce.all.filt.harmony@meta.data$celltype = factor(sce.all.filt.harmony@meta.data$celltype, 
                                             levels = c("NKTcell", "Tcell", "Bcell", "Monocyte/Macrophage", "ProNKcell", "NKcell"))

genes_to_check = c('CD3D', 'CD3E', 'CD4','CD8A', # T Cells 
                   'CD19', 'CD79A', 'MS4A1', # B cells
                   'IGHG1', 'MZB1', 'SDC1', # Plasma cells
                   'CD68', 'CD163', 'CD14', 'C1QA', 'C1QB', 'ITGAM', 'AIF1',# macrophages
                   'TPSAB1', 'TPSB2', # mast cells,
                   'CD14', 'S100A9', 'S100A8', 'MMP19', # monocyte
                   'FCGR3A', 'FGFBP2', 'CX3CR1', 'KLRB1', 'NCR1', 'NKG7', # NK cells
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'MKI67'
)


averageHeatmap(object = sce.all.filt.harmony, markerGene = genes_to_check, cluster_rows = T)

# 标记基因
Idents(sce.all.filt.harmony) = sce.all.filt.harmony$celltype
allCells = names(Idents(sce.all.filt.harmony))
allType = levels(Idents(sce.all.filt.harmony))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce.all.filt.harmony)== x ]
  cg = sample(cgCells, 500)
  cg
}))

cg_sce = sce.all.filt.harmony[, allCells %in% choose_Cells]
cg_sce
sce.markers = FindAllMarkers(object = cg_sce, min.pct = 0.25, only.pos = T)
sce.markers = FindAllMarkers(object = sce.all.filt.harmony, min.pct = 0, only.pos = T)

sce.markers$gene = rownames(sce.markers)
write.csv(sce.markers, file = "./3-celltype/1-allcell-marker.csv")

# new figures
library(scop)
CellDimPlot(
  srt = sce.all.filt.harmony, group.by = c("celltype", "Group"), palcolor = col_vec,
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave(filename = "./3-celltype/3-umap.celltype2.pdf", width = 10, height = 4.5)

FeatureDimPlot(
  srt = sce.all.filt.harmony,
  features = c("CD3D", "CD79A", "CD68", "FCGR3A", "NKG7", "MKI67"),
  reduction = "UMAP",
  theme_use = "theme_blank"
)
ggsave(filename = "./3-celltype/4-umap.markers.pdf", width = 10, height = 4.5)


features = c('ARHGAP17', 'BTK', 'NAGK', 'PAK1', 'PRR13', 'QPCT', 'RYBP', 'SLC22A4', 'TAGLN', 'TMEM154')
averageHeatmap(object = sce.all.filt.harmony, markerGene = features, cluster_rows = T)

FeatureDimPlot(
  srt = sce.all.filt.harmony,
  features = "QPCT",
  reduction = "UMAP",
  theme_use = "theme_blank"
)
ggsave(filename = "./3-celltype/4-umap.markers2.pdf", width = 5, height = 3.5)
