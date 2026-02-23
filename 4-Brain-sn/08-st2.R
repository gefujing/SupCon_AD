
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

load("./5-neusub/neusub2.Rdata")
load("./6-macrosub/microsub.Rdata")
load("./2-qc/1-sce.all.filt.Rdata")
load("./1-datapepare/mycol.Rdata")
load("./4-celltype/3-sce.all.filt.int.celltype.Rdata")

# 注释细胞
s100_norm = GetAssayData(sce.all.filt, slot = "counts")["S100A12", ]
s100_norm_new = s100_norm[colnames(microsub)]
summary(is.na(s100_norm_new))

microsub$celltype = ifelse(s100_norm_new > 0, "S100A12.pos", "S100A12.neg")
table(microsub$celltype)

S100A12_id = rownames(microsub@meta.data)[microsub$celltype == "S100A12.pos"]
sce.all.filt.int@meta.data$celltype = ifelse(rownames(sce.all.filt.int@meta.data) %in% S100A12_id,
                                             paste0("S100A12_", sce.all.filt.int@meta.data$celltype),
                                             sce.all.filt.int@meta.data$celltype)

damaged_id = rownames(neusub@meta.data)[neusub$damaged_neuron == "Damaged"]
sce.all.filt.int@meta.data$celltype = ifelse(rownames(sce.all.filt.int@meta.data) %in% damaged_id,
                                             paste0("Damaged_", sce.all.filt.int@meta.data$celltype),
                                             sce.all.filt.int@meta.data$celltype)

table(sce.all.filt.int@meta.data$celltype)

# 取子集
meta = sce.all.filt.int@meta.data
meta$cell = rownames(meta)

cells_keep <- meta %>%
  group_split(celltype) %>%               # 按 celltype 拆成 list(list of data.frames)
  lapply(function(df) {
    # 每个 cell type 内部抽样：如果少于100个，就全保留
    slice_sample(df, n = min(100, nrow(df)))
  }) %>%
  bind_rows() %>%                         # 再拼回一个 data.frame
  pull(cell)                              # 取出细胞名（barcode）


# 用这些细胞做子集
cg_sce = subset(sce.all.filt.int, cells = cells_keep)

# 
epimatrix = cg_sce@assays$RNA@layers$counts
rownames(epimatrix) = rownames(cg_sce)
colnames(epimatrix) = colnames(cg_sce)
write.csv(epimatrix, file = "./8-st/cellmatrix2.csv")

meta = cg_sce@meta.data[,c(1,18)]
write.csv(meta, file = "./8-st/cellmeta2.csv")

