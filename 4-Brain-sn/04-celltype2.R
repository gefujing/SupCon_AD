
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
load(file = "./4-celltype/3-sce.all.filt.int.celltype.Rdata")
load(file = "./1-datapepare/mycol.Rdata")

# 标记基因
Idents(sce.all.filt.int) = sce.all.filt.int@meta.data$celltype
sce.all.filt.int@meta.data$celltype = factor(sce.all.filt.int@meta.data$celltype, 
                                             levels = c("Astro", "Excit", "Inhibit", "Oligo", "OPC", "Endothelial", "VLMC", "Microglia"))

genes_to_check = c('AQP4', 'GJA1', 'ALDOC', 'GFAP', 'CLU', 'ALDH1L1', 'SLC1A3', 'GLUL', 'SLC1A2', # Astro
                   'PDGFRA', 'CSPG4', 'OLIG1', 'OLIG2', 'SOX10',  # OPC
                   'PLP1', 'MBP', 'MOBP', 'MOG', # Oligo
                   'SNAP25', 'STMN2', 'SYT1', 'SLC17A7', 'GRIN2A', 'SATB2', # Excit
                   'GAD1', 'GAD2', # Inhibit
                   'C1QA', 'C1QB', 'C1QC', 'TYROBP', 'P2RY12', 'HEXB', 'TREM2', 'CTSS', 'CSF1R', 'CD74',  # Mic
                   'CLDN5', 'NOSTRIN', 'FLT1', 'ITM2A', 'KLF2', 'BSG', # Endo
                   'DCN', 'COL6A2', 'ITIH5' # vascular & leptomeningeal cells
)


averageHeatmap(object = sce.all.filt.int, markerGene = genes_to_check, cluster_rows = T)

# 标记基因
Idents(sce.all.filt.int) = sce.all.filt.int$celltype
allCells = names(Idents(sce.all.filt.int))
allType = levels(Idents(sce.all.filt.int))

choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce.all.filt.int)== x ]
  cg = sample(cgCells, 500)
  cg
}))

cg_sce = sce.all.filt.int[, allCells %in% choose_Cells]
cg_sce
sce.markers = FindAllMarkers(object = cg_sce, min.pct = 0.25, only.pos = T)

sce.markers$gene = rownames(sce.markers)
write.csv(sce.markers, file = "./4-celltype/1-allcell-marker.csv")


# 糖基化特征
sce.all.filt.int1 = subset(sce.all.filt.int, gse_id == "GSE157827")
sce.all.filt.int2 = subset(sce.all.filt.int, gse_id == "GSE167494")
sce.all.filt.int3 = subset(sce.all.filt.int, gse_id == "GSE174367")

genes_to_check = c("PVR", "MAN2A2", "MGAT2", "MGAT3", "MGAT4A", "MGAT5", "B4GALT1",  "ST6GAL1", "FUT8")

DotPlot(sce.all.filt.int, features = unique(genes_to_check), group.by = "celltype", assay='RNA')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p1 = DotPlot(sce.all.filt.int, features = unique(genes_to_check), group.by = "celltype", split.by = "diagnosis", assay='RNA', cols = col_vec)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdata = p1@data
pdata$celltype = str_split(string = pdata$id, pattern = "_", simplify = T)[,1]
pdata$group = str_split(string = pdata$id, pattern = "_", simplify = T)[,2]
pdata$group = factor(pdata$group, levels = c("HC", "AD"))

p2 = ggplot(data = pdata, mapping = aes(x = group, y = features.plot, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  facet_wrap(.~celltype, ncol = 4)+
  theme_bw()+
  #  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))+
  scale_color_viridis()
p2
ggsave(p2, filename = "./4-celltype/10-Gly-genes.pdf", width = 6, height = 4.3)

features = c('ARHGAP17', 'BTK', 'NAGK', 'PAK1', 'PRR13', 'QPCT', 'RYBP', 'SLC22A4', 'TAGLN', 'TMEM154')
averageHeatmap(object = sce.all.filt.int, markerGene = features, cluster_rows = T)

ggsave(filename = "./3-celltype/4-umap.markers.pdf", width = 10, height = 4.5)

# new figures
CellDimPlot(
  srt = sce.all.filt.int, group.by = c("celltype"), split.by = "gse_id", palcolor = col_vec,
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave(filename = "./4-celltype/3-umap.celltype2.pdf", width = 12, height = 4.5)


