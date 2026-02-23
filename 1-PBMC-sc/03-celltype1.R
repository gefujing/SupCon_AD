
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

# 设置数据
load(file = "./2-harmony/sce.all.filt.harmony.Rdata")
load(file = "./1-qc/mycol.Rdata")

# 检测分群相关性
DimPlot(object = sce.all.filt.harmony, group.by = "seurat_clusters", reduction = "umap", label = T, cols = col_vec)
ggsave(filename = "./3-celltype/1-cluster.umap.pdf", width = 4.75, height = 3.75)

DimPlot(object = sce.all.filt.harmony, group.by = "orig.ident", reduction = "umap", cols = col_vec) 
DimPlot(object = sce.all.filt.harmony, group.by = "Group", reduction = "umap", cols = col_vec)

# 检查总体基因情况
genes_to_check = c('PTPRC', # immune cell
                   'CD3D', 'CD3E', 'CD4','CD8A', # T Cells 
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

genes_to_check
DotPlot(sce.all.filt.harmony, features = unique(genes_to_check), group.by = "seurat_clusters", assay='RNA') + coord_flip()
ggsave(filename = "./3-celltype/2-cluster.markers2.pdf", width = 8.5, height = 5.5)

# 细胞注释
celltype = data.frame(ClusterID = 0:25,
                      celltype = 0:25) 
## 定义细胞亚群
celltype[celltype$ClusterID %in% c(0,3,4,5,6,7,14,15,21,24), 2] = "Tcell"
celltype[celltype$ClusterID %in% c(2,9,7,25), 2] = "NKTcell"
celltype[celltype$ClusterID %in% c(1,11,13,17), 2] = "NKcell"
celltype[celltype$ClusterID %in% c(19), 2] = "ProNKcell"
celltype[celltype$ClusterID %in% c(8,12,18,23), 2] = "Bcell"
celltype[celltype$ClusterID %in% c(10,16,20,22), 2] = "Monocyte/Macrophage" 


## 写入细胞亚群
table(celltype$celltype)
sce.all.filt.harmony@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  sce.all.filt.harmony@meta.data[which(sce.all.filt.harmony@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]}
table(sce.all.filt.harmony@meta.data$celltype)

## 去除双细胞
Idents(sce.all.filt.harmony) = "seurat_clusters"
k = sce.all.filt.harmony@meta.data$seurat_clusters %in% c(17)
k = rownames(sce.all.filt.harmony@meta.data)[!k]

sce.all.filt.harmony = subset(sce.all.filt.harmony, cells = k)
save(sce.all.filt.harmony, file = "./3-celltype/sce.all.filt.harmony.celltype.Rdata")

# 查看细胞亚群
DimPlot(sce.all.filt.harmony, reduction = "umap", group.by = "celltype", label = T, cols = col_vec) 
ggsave(filename = "./3-celltype/3-umap.celltype.pdf", width = 6, height = 3.75)

DimPlot(sce.all.filt.harmony, reduction = "umap", group.by = "celltype", split.by = "Group", label = T, cols = col_vec) 
ggsave(filename = "./3-celltype/4-umap.celltype.group.pdf", width = 8, height = 3.75)



# 标记基因
Idents(sce.all.filt.harmony) = sce.all.filt.harmony$seurat_clusters
allCells = names(Idents(sce.all.filt.harmony))
allType = levels(Idents(sce.all.filt.harmony))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce.all.filt.harmony)== x ]
  cg = sample(cgCells, 50)
  cg
}))

cg_sce = sce.all.filt.harmony[, allCells %in% choose_Cells]
cg_sce
sce.markers = FindMarkers(object = cg_sce, ident.1 = 19, min.pct = 0.25, only.pos = T)
