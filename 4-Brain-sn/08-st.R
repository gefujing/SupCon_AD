
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
load("./4-celltype/3-sce.all.filt.int.celltype.Rdata")
load("./1-datapepare/mycol.Rdata")

# 注释细胞
epimatrix = microsub@assays$RNA@layers$counts
rownames(epimatrix) = rownames(microsub)
colnames(epimatrix) = colnames(microsub)

microsub$celltype = ifelse(epimatrix["QPCT",] > 0, "QPCT.pos", "QPCT.neg")
table(microsub$celltype)

QPCT_id = rownames(microsub@meta.data)[microsub$celltype == "QPCT.pos"]
sce.all.filt.int@meta.data$celltype = ifelse(rownames(sce.all.filt.int@meta.data) %in% QPCT_id,
                                             paste0("QPCT_", sce.all.filt.int@meta.data$celltype),
                                             sce.all.filt.int@meta.data$celltype)

damaged_id = rownames(neusub@meta.data)[neusub$damaged_neuron == "Damaged"]
sce.all.filt.int@meta.data$celltype = ifelse(rownames(sce.all.filt.int@meta.data) %in% damaged_id,
                                             paste0("Damaged_", sce.all.filt.int@meta.data$celltype),
                                             sce.all.filt.int@meta.data$celltype)

table(sce.all.filt.int@meta.data$celltype)

# 取子集
Idents(sce.all.filt.int) = sce.all.filt.int$celltype
allCells = names(Idents(sce.all.filt.int))
allType = levels(Idents(sce.all.filt.int))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce.all.filt.int)== x ]
  cg = sample(cgCells, 100)
  cg
}))

cg_sce = sce.all.filt.int[, allCells %in% choose_Cells]
cg_sce

# 
epimatrix = cg_sce@assays$RNA@layers$counts
rownames(epimatrix) = rownames(cg_sce)
colnames(epimatrix) = colnames(cg_sce)
write.csv(epimatrix, file = "./8-st/cellmatrix.csv")

meta = cg_sce@meta.data[,c(1,18)]
write.csv(meta, file = "./8-st/cellmeta.csv")
