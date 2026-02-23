
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

# 设置数据
load(file = "./3-harmony/2-sce.all.filt.int.Rdata")
load(file = "./1-datapepare/mycol.Rdata")

# 检测分群相关性
table(sce.all.filt.int$seurat_clusters)
av = AverageExpression(sce.all.filt.int, group.by = "seurat_clusters", assays = "RNA")
av = av[[1]]
cg = names(tail(sort(apply(av, 1, sd)), 2000))
kk = cor(as.matrix(av[cg,]), method = 'spearman')

ggheatmap(kk,
          scale = "none",
          cluster_rows = T,
          cluster_cols = T,
          border = "black",
          show_cluster_cols = F,
          show_cluster_rows = F,
          color = c("blue", "white", "red"))

ggsave(filename = "./4-celltype/1-cluster.cor.heatmap.pdf", width = 7.1, height = 6.1)

# 检查总体基因情况
genes_to_check = c('AQP4', 'GJA1', 'ALDOC', 'GFAP', 'CLU', 'ALDH1L1', 'SLC1A3', 'GLUL', 'SLC1A2', # Astro
                   'PDGFRA', 'CSPG4', 'OLIG1', 'OLIG2', 'SOX10',  # OPC
                   'PLP1', 'MBP', 'MOBP', 'MOG', # Oligo
                   'SNAP25', 'STMN2', 'SYT1', 'SLC17A7', 'GRIN2A', 'SATB2', # Excit
                   'GAD1', 'GAD2', # Inhibit
                   'C1QA', 'C1QB', 'C1QC', 'TYROBP', 'P2RY12', 'HEXB', 'TREM2', 'CTSS', 'CSF1R', 'CD74',  # Mic
                   'CLDN5', 'NOSTRIN', 'FLT1', 'ITM2A', 'KLF2', 'BSG', # Endo
                   'COLA1A2', 'DCN', 'LUM', 'COL6A2', 'ITIH5' # vascular & leptomeningeal cells
)

genes_to_check
p <- DotPlot(sce.all.filt.int, features = unique(genes_to_check), group.by = "seurat_clusters", assay='RNA') + coord_flip()
p
ggsave(filename = "./4-celltype/2-cluster.markers.pdf", width = 13, height = 8)

# 细胞注释
celltype = data.frame(ClusterID = 0:43, celltype = 0:43) 

## 定义细胞亚群
celltype[celltype$ClusterID %in% c(1,4,7,11,23,36), 2] = 'Astro'
celltype[celltype$ClusterID %in% c(3), 2] = 'OPC'
celltype[celltype$ClusterID %in% c(0,6,10,12,14,29,33,34,41), 2] = 'Oligo'
celltype[celltype$ClusterID %in% c(2,5,9,15,17,19,24,26,27,31,32,35,37,38,39,40,42), 2] = 'Excit'
celltype[celltype$ClusterID %in% c(13,16,18,20,22,28), 2] = 'Inhibit'
celltype[celltype$ClusterID %in% c(8,21), 2] = 'Microglia'
celltype[celltype$ClusterID %in% c(30), 2] = 'Endothelial' 
celltype[celltype$ClusterID %in% c(25,43), 2] = 'VLMC'


## 写入细胞亚群
table(celltype$celltype)
sce.all.filt.int@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  sce.all.filt.int@meta.data[which(sce.all.filt.int@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]}

table(sce.all.filt.int@meta.data$celltype)
save(sce.all.filt.int, file = "./4-celltype/3-sce.all.filt.int.celltype.Rdata")

# 查看细胞亚群
p1 = DimPlot(sce.all.filt.int, reduction = "umap", group.by = "celltype", cols = col_vec, label = T) + NoLegend()
p2 = DimPlot(sce.all.filt.int, reduction = "umap", group.by = "celltype", split.by = "gse_id", cols = col_vec)  + NoLegend()
pall = p1 + p2 + plot_layout(ncol = 2, widths = c(1, 3))
ggsave(filename = "./4-celltype/3-umap-celltype.pdf", width = 12.5, height = 4)




# 标记基因
Idents(sce.all.filt.int) = sce.all.filt.int$seurat_clusters
allCells = names(Idents(sce.all.filt.int))
allType = levels(Idents(sce.all.filt.int))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce.all.filt.int)== x ]
  cg = sample(cgCells, 100)
  cg
}))

cg_sce = sce.all.filt.int[, allCells %in% choose_Cells]
cg_sce
sce.markers = FindMarkers(object = cg_sce, ident.1 = 24, min.pct = 0.25, only.pos = T)
sce.markers$gene = rownames(sce.markers)


