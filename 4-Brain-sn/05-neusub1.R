

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
neusub = subset(x = sce.all.filt.int, celltype %in% c("Excit", "Inhibit"))
table(neusub@meta.data$orig.ident)

# 取子集
# neusub = subset(x = neusub, seurat_clusters %in% c(0:8, 10:17, 19:20, 22:29))

## 创建Seurat数据
counts = neusub@assays$RNA@layers$counts
rownames(counts) = rownames(neusub)
colnames(counts) = colnames(neusub)

meta.data = neusub@meta.data
neusub = CreateSeuratObject(counts = counts, meta.data = meta.data)
rm(list = c("counts", "meta.data", "sce.all.filt.int"))

## 添加分组信息
table(neusub@meta.data$orig.ident)
table(neusub@meta.data$diagnosis)

## 标准化，归一化，PCA
neusub = neusub[!grepl("^RP[SL]", rownames(neusub), ignore.case = T), ]
neusub = NormalizeData(neusub, normalization.method =  "LogNormalize", scale.factor = 1e4)
neusub = FindVariableFeatures(neusub, selection.method = "vst", nfeatures = 2000) 
neusub = ScaleData(neusub) 
neusub = RunPCA(object = neusub, features = VariableFeatures(neusub))
neusub = RunHarmony(neusub, group.by.vars = "orig.ident")
neusub = RunUMAP(neusub, dims = 1:20, reduction = "harmony")

neusub = FindNeighbors(neusub, reduction = "harmony", dims = 1:20) 
neusub = FindClusters(neusub, resolution = 0.5, algorithm = 1)
table(neusub@meta.data$seurat_clusters)

DimPlot(neusub, reduction = "umap", group.by = "seurat_clusters", label = T, cols = col_vec)
ggsave(filename = "./5-neusub/1-neusub.umpa.pdf", width = 5.5, height = 4) 


# celltype
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
p <- DotPlot(neusub, features = unique(genes_to_check), group.by = "seurat_clusters", assay='RNA') + coord_flip() + theme_bw()
p
ggsave(filename = "./5-neusub/2-cluster.markers.pdf", width = 8.5, height = 6.5)


# 细胞注释
celltype = data.frame(ClusterID = 0:26,
                      celltype = 0:26) 
## 定义细胞亚群
celltype[celltype$ClusterID %in% c(0,2,3,4,5,8,9,11,12,15,16,17,18,19,20,22,24,25), 2] = "Excit"
celltype[celltype$ClusterID %in% c(1,6,7,10,13,14,21,23,26), 2] = "Inhibit"

## 写入细胞亚群
table(celltype$celltype)
neusub@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  neusub@meta.data[which(neusub@meta.data$seurat_clusters == celltype$ClusterID[i]), 'celltype'] <- celltype$celltype[i]}
table(neusub@meta.data$celltype)

neusub@meta.data$celltype = factor(neusub@meta.data$celltype,
                                   levels = c("Excit", "Inhibit"))
save(neusub, file = "./5-neusub/neusub.Rdata")


DimPlot(neusub, reduction = "umap", group.by = "celltype", label = T, cols = col_vec)
ggsave(filename = "./5-neusub/1-neusub.umpa.cell.pdf", width = 4.5, height = 3.5) 


# AUCell
exprMatrix = neusub@assays$RNA@layers$counts
rownames(exprMatrix) = rownames(neusub)
colnames(exprMatrix) = colnames(neusub)

geneset = list(
  damaged = c("ATF3",
              "JUN","SOX11","STAT3","KLF6",
              # "GAP43","SPRR1A","CAP23","ECEL1",
              # "BAX","BCL2","CASP3","HSP27",
              # "VAMP4","SYT11","NEFH","MAP2",
              # "IL6","TNF","CXCL10","GFAP","TREM2","C1QA",
              "MICA","MICB","ULBP1","ULBP2","ULBP3","ULBP4","ULBP5","ULBP6")
)

cells_AUC = AUCell_run(exprMatrix, geneset)

damaged_score = as.numeric(getAUC(cells_AUC)["damaged", ])
neusub@meta.data$damaged_score = damaged_score

DotPlot(neusub, features = c("damaged_score"), group.by = "diagnosis")
ggsave(filename = "./5-neusub/3-damage.pdf", width = 3, height = 1.4)
FeaturePlot(neusub, features = c("damaged_score"), split.by = "diagnosis")


neusub$damaged_neuron = ifelse(neusub$damaged_score > median(neusub$damaged_score), "Damaged", "Normal")
neusub$damaged_neuron = factor(neusub$damaged_neuron, levels = c("Normal", "Damaged"))

neusub = RunDEtest(srt = neusub, group_by = "damaged_neuron", fc.threshold = 1, only.pos = FALSE)

VolcanoPlot(srt = neusub, group_by = "damaged_neuron")
ggsave(filename = "./5-neusub/4-damage2.pdf", width = 7, height = 3.8)

degs = neusub@tools[["DEtest_damaged_neuron"]][["AllMarkers_wilcox"]]
degs = degs[degs$group1 == "Damaged", ]

cg = degs$gene[(degs$p_val_adj < 1e-5) & (degs$avg_log2FC > 0.5)]
mm = read.csv(file = "./5-neusub/S2_File.csv")
mm = mm$ENTREZ.gene.symbol

gene = intersect(cg, mm)


load(file = "./4-celltype/3-sce.all.filt.int.celltype.Rdata")


damaged_id = rownames(neusub@meta.data)[neusub$damaged_neuron == "Damaged"]

sce.all.filt.int@meta.data$celltype = ifelse(rownames(sce.all.filt.int@meta.data) %in% damaged_id,
                                             paste0("Damaged_", sce.all.filt.int@meta.data$celltype),
                                             sce.all.filt.int@meta.data$celltype)

gene = c("MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "ULBP4", "ULBP5", "ULBP6")

DotPlot(neusub, features = gene, group.by = "damaged_neuron")+ coord_flip()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


ggsave(filename = "./5-neusub/4-damage4.pdf", width = 3.6, height = 2)







# 标记基因
Idents(neusub) = neusub$celltype
allCells = names(Idents(neusub))
allType = levels(Idents(neusub))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(neusub)== x ]
  cg = sample(cgCells, 500)
  cg
}))

cg_sce = neusub[, allCells %in% choose_Cells]
cg_sce
sce.markers = FindAllMarkers(object = neusub, min.pct = 0.25, only.pos = T)
write.csv(sce.markers, file = "./5-neusub/1-allcell-marker.csv")


# 可视化
# new figures
library(scop)
CellDimPlot(
  srt = neusub, group.by = c("celltype", "Group"), palcolor = col_vec,
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave(filename = "./5-neusub/3-umap.celltype2.pdf", width = 8, height = 4)

FeatureDimPlot(
  srt = neusub,
  features = c("CD14", "LYZ", "S100A12", "LGALS2", "IL1B", "HLA-DPB1", "MX1", "FCGR3A", "CD1C", "IRF7"),
  reduction = "UMAP",
  theme_use = "theme_blank",
  ncol = 5
)
ggsave(filename = "./5-neusub/4-umap.markers.pdf", width = 11, height = 5)



ht <- GroupHeatmap(
  srt = neusub,
  features = c("CD14", "LYZ", "S100A12", "LGALS2", "IL1B", "HLA-DPB1", "MX1", "FCGR3A", "CD1C", "IRF7"),
  group.by = c("celltype"),
  heatmap_palette = "YlOrRd",
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  row_names_side = "left",
  add_dot = TRUE)
print(ht$plot)









