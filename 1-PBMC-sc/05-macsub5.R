
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(tidyverse)
library(scop)
library(Seurat)
library(patchwork)

# 设置数据
load(file = "./5-macsub/macsub.Rdata")
load(file = "./1-qc/mycol.Rdata")
table(macsub$celltype)

# 差异基因分析
p_adj_cutoff = 0.2
min_pct = 0.05
logfc_thresh = 0
test_use = "wilcox"
min_cells_per_group = 1
assay_to_use = DefaultAssay(macsub)

obj = macsub

# 基本检查
stopifnot("celltype"  %in% colnames(obj@meta.data))
stopifnot("Group" %in% colnames(obj@meta.data))

DefaultAssay(obj) = assay_to_use
celltypes = sort(unique(obj@meta.data$celltype))

# 结果容器
summary_list = list()
markers_all = list()

for (ct in celltypes) {
  message("Processing cell type: ", ct)
  sub = subset(obj, subset = celltype == ct)
  diag_tab = table(sub@meta.data$Group)
  
  Idents(sub) = "Group"
  
  mk = tryCatch(
    FindMarkers(
      sub,
      ident.1 = "AD",
      ident.2 = "HC",
      test.use = test_use,
      min.pct = min_pct,
      logfc.threshold = logfc_thresh,
      only.pos = FALSE)
  )
  
  
  if (!("avg_log2FC" %in% colnames(mk))) {
    if ("avg_logFC" %in% colnames(mk)) {
      mk$avg_log2FC = mk$avg_logFC
    } else {
      stop("No avg_log2FC/avg_logFC columns...")
    }
  }
  
  mk = mk %>%
    rownames_to_column(var = "gene") %>%
    mutate(signif = p_val_adj < p_adj_cutoff,
           direction = case_when(
             signif & avg_log2FC > 0  ~ "up_in_AD",
             signif & avg_log2FC < 0  ~ "down_in_AD",
             TRUE                     ~ "ns"
           ))
  
  n_up   = sum(mk$direction == "up_in_AD",   na.rm = TRUE)
  n_down = sum(mk$direction == "down_in_AD", na.rm = TRUE)
  n_sig  = n_up + n_down
  
  summary_list[[ct]] = tibble(
    celltype   = ct,
    n_up_in_AD = n_up,
    n_down_in_AD = n_down,
    n_sig_total = n_sig,
    n_cells_AD = as.integer(diag_tab["AD"]),
    n_cells_HC = as.integer(diag_tab["HC"])
  )
  
  markers_all[[ct]] = mk
}


for (ct in names(markers_all)) {
  fn <- paste0("./5-macsub/14-DEGs/3-DEGs_", gsub("[^A-Za-z0-9_]+","_", ct), "_AD_vs_HC.csv")
  write.csv(markers_all[[ct]], fn, row.names = FALSE)
}
