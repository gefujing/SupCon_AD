
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
library(AUCell)
load("./4-celltype/3-sce.all.filt.int.celltype.Rdata")
load(file = "./1-datapepare/mycol.Rdata")

# 亚群再分析
table(sce.all.filt.int@meta.data$celltype)
microsub = subset(x = sce.all.filt.int, celltype %in% c("Microglia"))
table(microsub@meta.data$orig.ident)

## 创建Seurat数据
counts = microsub@assays$RNA@layers$counts
rownames(counts) = rownames(microsub)
colnames(counts) = colnames(microsub)

meta.data = microsub@meta.data
microsub = CreateSeuratObject(counts = counts, meta.data = meta.data)
rm(list = c("counts", "meta.data", "sce.all.filt.int"))

## 添加分组信息
table(microsub@meta.data$orig.ident)
table(microsub@meta.data$diagnosis)

## 标准化，归一化，PCA
microsub = microsub[!grepl("^RP[SL]", rownames(microsub), ignore.case = T), ]
microsub = NormalizeData(microsub, normalization.method =  "LogNormalize", scale.factor = 1e4)
microsub = FindVariableFeatures(microsub, selection.method = "vst", nfeatures = 2000) 
microsub = ScaleData(microsub) 
microsub = RunPCA(object = microsub, features = VariableFeatures(microsub))
microsub = RunHarmony(microsub, group.by.vars = "orig.ident")
microsub = RunUMAP(microsub, dims = 1:20, reduction = "harmony")

microsub = FindNeighbors(microsub, reduction = "harmony", dims = 1:20) 
microsub = FindClusters(microsub, resolution = 0.5, algorithm = 1)
table(microsub@meta.data$seurat_clusters)

microsub = subset(x = microsub, seurat_clusters %in% c(0:8, 10:12))

DimPlot(microsub, reduction = "umap", group.by = "seurat_clusters", label = T, cols = col_vec)
ggsave(filename = "./6-microsub/1-microsub.umpa.pdf", width = 5.5, height = 4) 

save(microsub, file = "./6-microsub/microsub.Rdata")




# 标记基因
genes_to_check = c('C1QA', 'C1QB', 'C1QC', 'TYROBP', 'P2RY12', 'HEXB', 'TREM2', 'CTSS', 'CSF1R', 'CD74')

genes_to_check
p <- DotPlot(microsub, features = unique(genes_to_check), group.by = "seurat_clusters", assay='RNA') + coord_flip()
p
ggsave(filename = "./4-celltype/2-cluster.markers.pdf", width = 13, height = 8)

FeatureDimPlot(srt = microsub, features = genes_to_check, ncol = 5)




