

# 准备环境
rm(list=ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(stringr)
library(clustree)
library(cowplot)
library(dplyr)
library(stringr)
library(data.table)
library(patchwork)
library(ggsci)
library(harmony)
library(scRNAtoolVis)
library(scop)
load("./3-celltype/sce.all.filt.harmony.celltype.Rdata")
load(file = "./1-qc/mycol.Rdata")

# 亚群再分析
table(macsub@meta.data$celltype)
macsub = subset(x = macsub, celltype == "Monocyte/Macrophage")
table(macsub@meta.data$orig.ident)

# 取子集
# macsub = subset(x = macsub, seurat_clusters %in% c(0:6, 14))

## 创建Seurat数据
counts = macsub@assays$RNA@layers$counts
rownames(counts) = rownames(macsub)
colnames(counts) = colnames(macsub)

meta.data = macsub@meta.data
macsub = CreateSeuratObject(counts = counts, meta.data = meta.data)
rm(list = c("counts", "meta.data", "macsub"))

## 添加分组信息
table(macsub@meta.data$orig.ident)
table(macsub@meta.data$Group)

## 标准化，归一化，PCA
macsub = macsub[!grepl("^RP[SL]", rownames(macsub), ignore.case = T), ]
macsub = NormalizeData(macsub, normalization.method =  "LogNormalize", scale.factor = 1e4)
macsub = FindVariableFeatures(macsub, selection.method = "vst", nfeatures = 2000) 
macsub = ScaleData(macsub) 
macsub = RunPCA(object = macsub, features = VariableFeatures(macsub))
macsub = RunHarmony(macsub, group.by.vars = "orig.ident")
macsub = RunUMAP(macsub, dims = 1:20, reduction = "harmony")

macsub = FindNeighbors(macsub, reduction = "harmony", dims = 1:20) 
macsub = FindClusters(macsub, resolution = 0.5, algorithm = 1)
table(macsub@meta.data$seurat_clusters)

DimPlot(macsub, reduction = "umap", group.by = "seurat_clusters", label = T, cols = col_vec)
ggsave(filename = "./5-macsub/1-macsub.umpa.pdf", width = 4.2, height = 4) 


# celltype
genes_to_check = c('CD14', 'LYZ', 'S100A8','S100A9', 'S100A12', # CD14 Mono 
                   'FCGR3A', 'CX3CR1', 'MS4A7', 'LILRB2', 'IFITM3', # CD16 Mono
                   'CLEC10A', 'HLA-DRA', 'BATF3', 'IRF8', 'CADM1', # cDC1
                   'CD1C', 'FCER1A', 'CLEC10A', 'ITGAX', # cDC2
                   'CLEC4C', 'GZMB', 'JCHAIN', 'TCF4', 'IRF7', # pDC
                   'CEACAM8', 'FCGR3B', 'CSF3R', 'MPO', 'CEACAM8', # Neutrophils
                   'C1QA', 'C1QB', 'CD68', 'CD163' # Macrophage
)


genes_to_check
p <- DotPlot(macsub, features = unique(genes_to_check), group.by = "seurat_clusters", assay='RNA') + coord_flip() + theme_bw()
p
ggsave(filename = "./5-macsub/2-cluster.markers.pdf", width = 4.5, height = 5)


# 细胞注释
celltype = data.frame(ClusterID = 0:7,
                      celltype = 0:7) 
## 定义细胞亚群
celltype[celltype$ClusterID %in% c(1), 2] = "S100A12+CD14+Mono"
celltype[celltype$ClusterID %in% c(2), 2] = "LGALS2+CD14+Mono"
celltype[celltype$ClusterID %in% c(3), 2] = "IL1B+CD14+Mono"
celltype[celltype$ClusterID %in% c(4), 2] = "HLA-DPB1+CD14+Mono"
celltype[celltype$ClusterID %in% c(6), 2] = "MX1+CD14+Mono"
celltype[celltype$ClusterID %in% c(0), 2] = "CD16+Mono" 
celltype[celltype$ClusterID %in% c(5), 2] = "cDC2" 
celltype[celltype$ClusterID %in% c(7), 2] = "pDC" 

## 写入细胞亚群
table(celltype$celltype)
macsub@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  macsub@meta.data[which(macsub@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]}
table(macsub@meta.data$celltype)

macsub@meta.data$celltype = factor(macsub@meta.data$celltype,
                                   levels = c("S100A12+CD14+Mono", "LGALS2+CD14+Mono", "IL1B+CD14+Mono", "HLA-DPB1+CD14+Mono",
                                              "MX1+CD14+Mono", "CD16+Mono", "cDC2", "pDC"))
save(macsub, file = "./5-macsub/macsub.Rdata")


DimPlot(macsub, reduction = "umap", group.by = "celltype", label = T, cols = col_vec)
ggsave(filename = "./5-macsub/1-macsub.umpa.cell.pdf", width = 4.5, height = 3.5) 


# 标记基因
Idents(macsub) = macsub$celltype
allCells = names(Idents(macsub))
allType = levels(Idents(macsub))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(macsub)== x ]
  cg = sample(cgCells, 500)
  cg
}))

cg_sce = macsub[, allCells %in% choose_Cells]
cg_sce
sce.markers = FindAllMarkers(object = macsub, min.pct = 0.25, only.pos = T)
write.csv(sce.markers, file = "./5-macsub/1-allcell-marker.csv")


# 可视化
# new figures
library(scop)
CellDimPlot(
  srt = macsub, group.by = c("celltype", "Group"), palcolor = col_vec,
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave(filename = "./5-macsub/3-umap.celltype2.pdf", width = 8, height = 4)

FeatureDimPlot(
  srt = macsub,
  features = c("CD14", "LYZ", "S100A12", "LGALS2", "IL1B", "HLA-DPB1", "MX1", "FCGR3A", "CD1C", "IRF7"),
  reduction = "UMAP",
  theme_use = "theme_blank",
  ncol = 5
)
ggsave(filename = "./5-macsub/4-umap.markers.pdf", width = 11, height = 5)



ht <- GroupHeatmap(
  srt = macsub,
  features = c("CD14", "LYZ", "S100A12", "LGALS2", "IL1B", "HLA-DPB1", "MX1", "FCGR3A", "CD1C", "IRF7"),
  group.by = c("celltype"),
  heatmap_palette = "YlOrRd",
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  row_names_side = "left",
  add_dot = TRUE)
print(ht$plot)


features = c('ARHGAP17', 'BTK', 'NAGK', 'PAK1', 'PRR13', 'QPCT', 'RYBP', 'SLC22A4', 'TAGLN', 'TMEM154')

averageHeatmap(object = macsub, markerGene = features, cluster_rows = T, cluster_columns = T)






