
# 准备环境
rm(list = ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
library(harmony)
load(file = "./2-qc/1-sce.all.filt.Rdata")
load(file = "./1-datapepare/mycol.Rdata")

# 取子集
table(sce.all.filt.int@meta.data$seurat_clusters)
sce.all.filt = subset(x = sce.all.filt.int, seurat_clusters %in% c(0:30,32:34,36,38:39,41:43,45:49))
rm(list = c("sce.all.filt.int"))

# 标准化数据
sce.all.filt = NormalizeData(sce.all.filt, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce.all.filt = FindVariableFeatures(sce.all.filt)
sce.all.filt = ScaleData(sce.all.filt)
sce.all.filt = RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))
sce.all.filt = RunHarmony(sce.all.filt, group.by.vars = "orig.ident")
sce.all.filt = RunUMAP(sce.all.filt, dims = 1:20, reduction = "harmony")

DimPlot(object = sce.all.filt, group.by = "orig.ident", split.by = "gse_id", reduction = "umap", cols = col_vec)
ggsave(filename = "./3-harmony/1-dim-orig-spilt.pdf", width = 15, height = 4)

# 细胞聚类
sce.all.filt = FindNeighbors(sce.all.filt, reduction = "harmony", dims = 1:20)
sce.all.filt = FindClusters(sce.all.filt, resolution = 1, algorithm = 1)

DimPlot(sce.all.filt, reduction = "umap", group.by = "orig.ident", cols = col_vec, label = F)
ggsave(filename = "./3-harmony/1-dim-orig.pdf", width = 9.5, height = 3.75)

DimPlot(sce.all.filt, reduction = "umap", group.by = "seurat_clusters", cols = col_vec, label = F)
ggsave(filename = "./3-harmony/2-dim-cluster.pdf", width = 5.4, height = 3.75)

DimPlot(sce.all.filt, reduction = "umap", group.by = "seurat_clusters", cols = col_vec, label = T)
ggsave(filename = "./3-harmony/2-dim-cluster2.pdf", width = 5.4, height = 3.75)

# 保存数据
sce.all.filt = SetIdent(sce.all.filt, value = "seurat_clusters")
table(sce.all.filt@active.ident)
sce.all.filt.int = sce.all.filt
save(sce.all.filt.int, file = "./3-harmony/2-sce.all.filt.int.Rdata")
