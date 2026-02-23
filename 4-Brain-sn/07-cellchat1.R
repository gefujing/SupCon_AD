
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
library(ggheatmap)
library(ggsci)
library(org.Hs.eg.db)
library(CellChat)
library(scRNAtoolVis)

# 设置数据
load(file = "./4-scRNA/3-sce.all.filt.int.celltype.Rdata")
load(file = "./0-RawData/my_color.Rdata")

# 分组
Idents(sce.all.filt.int) = sce.all.filt.int@meta.data$celltype
sce.all.filt.int@meta.data$celltype = as.character(sce.all.filt.int@meta.data$celltype)

# episub
episub = sce.all.filt.int[ ,Idents(sce.all.filt.int) %in% c("PCC")]
table(sce.all.filt.int$celltype)

epimatrix = episub@assays$RNA@layers$counts
rownames(epimatrix) = rownames(episub)
colnames(epimatrix) = colnames(episub)

episub$celltype = ifelse(epimatrix["IGF2BP2",] > mean(epimatrix["IGF2BP2",]), "IGF2BP2.pos.PCC", "IGF2BP2.neg.PCC")
# episub$celltype = ifelse(epimatrix["IGF2BP2",] > 2, "IGF2BP2.pos.PCC", "IGF2BP2.neg.PCC")
table(episub$celltype)

sce.all.filt.int$celltype[match(colnames(episub), colnames(sce.all.filt.int))] = episub$celltype 
table(sce.all.filt.int$celltype)

sce.all.filt.int@meta.data$celltype = factor(sce.all.filt.int@meta.data$celltype, levels = c('IGF2BP2.pos.PCC', 'IGF2BP2.neg.PCC', 'Acinar', 'Fibroblast', 'Stellate', 'Endothelial', 'Macrophage', 'DCs', 'Tcell', 'Bcell', 'Plasma'))
DimPlot(sce.all.filt.int, reduction = "tsne", group.by = 'celltype', label = TRUE, repel = T, pt.size = 1, cols = col_vector)
save(episub, file = "./4-scRNA/episub.IGF2BP2.Rdata")


# fibsub
fibsub = sce.all.filt.int[ ,Idents(sce.all.filt.int) %in% c("Fibroblast")]
table(sce.all.filt.int$celltype)

fibmatrix = fibsub@assays$RNA@layers$counts
rownames(fibmatrix) = rownames(fibsub)
colnames(fibmatrix) = colnames(fibsub)

fibsub$celltype = ifelse(fibmatrix["ACTA2",] > mean(fibmatrix["ACTA2",]), "ACTA2.pos.FIB", "ACTA2.neg.FIB")
table(fibsub$celltype)

sce.all.filt.int@meta.data$celltype = as.character(sce.all.filt.int@meta.data$celltype)
sce.all.filt.int$celltype[match(colnames(fibsub), colnames(sce.all.filt.int))] = fibsub$celltype 
table(sce.all.filt.int$celltype)

sce.all.filt.int@meta.data$celltype = factor(sce.all.filt.int@meta.data$celltype, levels = c('IGF2BP2.pos.PCC', 'IGF2BP2.neg.PCC', "ACTA2.pos.FIB", "ACTA2.neg.FIB", 'Acinar', 'Stellate', 'Endothelial', 'Macrophage', 'DCs', 'Tcell', 'Bcell', 'Plasma'))
DimPlot(sce.all.filt.int, reduction = "tsne", group.by = 'celltype', label = TRUE, repel = T, pt.size = 1, cols = col_vector)
save(fibsub, file = "./4-scRNA/fibsub.ACTA2.Rdata")


# save data
save(sce.all.filt.int, file = "./4-scRNA/4-sce.all.filt.int.IGF2BP2.Rdata")


# write csv
cellmatrix = as.data.frame(sce.all.filt.int@assays$RNA@layers$counts)
rownames(cellmatrix) = rownames(sce.all.filt.int)
colnames(cellmatrix) = colnames(sce.all.filt.int)

write.csv(cellmatrix, gzfile("./4-scRNA/cellmatrix.csv.gz"), row.names = TRUE)


metadata = as.data.frame(sce.all.filt.int@meta.data)
write.csv(metadata, gzfile("./4-scRNA/cellmeta.csv.gz"), row.names = TRUE)








