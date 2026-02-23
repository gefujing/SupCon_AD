
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
library(pheatmap)
library(ggheatmap)
library(RColorBrewer)
library(ggthemes)
library(patchwork)
library(scRNAtoolVis)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(ggradar)

# 设置数据
load(file = "./4-celltype/3-sce.all.filt.int.celltype.Rdata")
load(file = "./1-datapepare/mycol.Rdata")
table(sce.all.filt.int$celltype)

# 差异基因分析
p_adj_cutoff = 0.05
min_pct = 0.10
logfc_thresh = 0
test_use = "wilcox"
min_cells_per_group = 1
assay_to_use = DefaultAssay(sce.all.filt.int)

obj = sce.all.filt.int

# 基本检查
stopifnot("celltype"  %in% colnames(obj@meta.data))
stopifnot("diagnosis" %in% colnames(obj@meta.data))

DefaultAssay(obj) = assay_to_use
celltypes = sort(unique(obj@meta.data$celltype))

# 结果容器
summary_list = list()
markers_all = list()

for (ct in celltypes) {
  message("Processing cell type: ", ct)
  sub = subset(obj, subset = celltype == ct)
  diag_tab = table(sub@meta.data$diagnosis)
  
  Idents(sub) = "diagnosis"
  
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

# 合并汇总
deg_summary <- bind_rows(summary_list) %>%
  arrange(desc(n_sig_total))

# 查看与保存
print(deg_summary)

# 可选：写出结果
write.csv(deg_summary, "./4-celltype/2-DEG_summary_by_celltype_AD_vs_HC.csv", row.names = FALSE)

for (ct in names(markers_all)) {
  fn <- paste0("./4-celltype/3-DEGs_", gsub("[^A-Za-z0-9_]+","_", ct), "_AD_vs_HC.csv")
  write.csv(markers_all[[ct]], fn, row.names = FALSE)
}

# 可视化
# 整理作图数据：上调为正值、下调为负值
plot_df = deg_summary %>%
  mutate(n_up_in_AD   = as.integer(n_up_in_AD),
         n_down_in_AD = as.integer(n_down_in_AD),
         n_sig_total  = as.integer(n_sig_total)) %>%
  transmute(celltype,
            n_sig_total,
            up   = n_up_in_AD,
            down = -n_down_in_AD  ) %>%
  pivot_longer(cols = c(up, down), names_to = "direction", values_to = "count_signed") %>%
  mutate(direction = recode(direction, up = "Up in AD", down = "Down in AD")) %>%
  mutate(celltype = fct_reorder(celltype, n_sig_total, .desc = TRUE))

plot_df_lab <- dplyr::filter(plot_df, count_signed != 0)

ggplot(plot_df, aes(x = celltype, y = count_signed, fill = direction)) +
  geom_col(width = 0.7) +
  geom_text(data = plot_df_lab,  aes(label = abs(count_signed)),  position = position_stack(vjust = 0.5),  size = 3.3) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  labs(x = "Cell type", y = "Number of DEGs (AD vs HC)") +
  theme_minimal(base_size = 12) +
  theme_classic() +
  scale_fill_manual(values = col_vec)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(filename = "./4-celltype/7-degs-summary.pdf", width = 5.5, height = 4)


# 基因功能富集分析
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggrepel)
library(ggtext)
library(fmsb)
library(scales)
library(stringr)
library(forcats)
library(glue)

# 检查数据
stopifnot(is.list(markers_all))
celltypes <- names(markers_all)
outdir <- "./4-celltype/8-enrich_outputs_go"
dir.create(outdir, showWarnings = FALSE)

# 设置参数
GO_ONT <- c("BP")      # 只做生物过程；想一起做：c("BP","CC","MF")
fdr_cutoff <- 0.05     # DE基因阈值
min_sig_genes <- 10    # 每组至少多少DEG才做富集
minGS <- 0; maxGS <- 1000
top_k_show <- 10       # 每个(细胞×方向×本体)图里展示多少条term
dir_cols <- c("Up in AD"="#4B4EAD","Down in AD"="#B31237")

## ===== 工具函数 =====
map_symbols_to_entrez <- function(genes){
  genes <- unique(na.omit(genes))
  tt <- suppressMessages(bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db))
  setNames(tt$ENTREZID, tt$SYMBOL)
}

do_go_ora <- function(gene_symbols, ct, dir_lab, ont){
  if(length(gene_symbols) < 5) return(NULL)
  g2e <- map_symbols_to_entrez(gene_symbols); if(length(g2e) < 5) return(NULL)
  ego <- tryCatch({
    enrichGO(gene = unname(g2e), OrgDb=org.Hs.eg.db, keyType="ENTREZID",
             ont = ont, pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.05,
             minGSSize=minGS, maxGSSize=maxGS, readable=TRUE)
  }, error=function(e) NULL)
  if(is.null(ego) || nrow(as.data.frame(ego))==0) return(NULL)
  as.data.frame(ego) |>
    transmute(source = ont, term = Description, ID,
              p.adjust, Count, GeneRatio, gene_symbols = geneID) |>
    mutate(celltype = ct, direction = dir_lab, score = -log10(p.adjust)) |>
    arrange(p.adjust)
}

# 循环富集
enrich_all <- list()

for(ct in celltypes){
  df <- markers_all[[ct]]
  stopifnot(all(c("gene","p_val_adj","direction") %in% colnames(df)))
  up_genes   <- df |> filter(direction=="up_in_AD") |> pull(gene) |> unique()
  down_genes <- df |> filter(direction=="down_in_AD") |> pull(gene) |> unique()
  
  for(ont in GO_ONT){
    if(length(up_genes)   >= min_sig_genes){
      up_tbl <- do_go_ora(up_genes, ct, "Up in AD", ont)
      if(!is.null(up_tbl)) enrich_all[[length(enrich_all)+1]] <- up_tbl
    }
    if(length(down_genes) >= min_sig_genes){
      dn_tbl <- do_go_ora(down_genes, ct, "Down in AD", ont)
      if(!is.null(dn_tbl)) enrich_all[[length(enrich_all)+1]] <- dn_tbl
    }
  }
}

enrich_df <- dplyr::bind_rows(enrich_all)
write.csv(enrich_df, file=file.path(outdir, "1-GO_enrichment_all.csv"), row.names=FALSE)


# 可视化
Excit_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Excit"), ]
Excit_up = Excit_up[!duplicated(Excit_up$ID),]; rownames(Excit_up) = Excit_up$ID

Astro_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Astro"), ]
Astro_up = Astro_up[!duplicated(Astro_up$ID),]; rownames(Astro_up) = Astro_up$ID

Oligo_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Oligo"), ]
Oligo_up = Oligo_up[!duplicated(Oligo_up$ID),]; rownames(Oligo_up) = Oligo_up$ID

Microglia_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Microglia"), ]
Microglia_up = Microglia_up[!duplicated(Microglia_up$ID),]; rownames(Microglia_up) = Microglia_up$ID

Inhibit_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Inhibit"), ]
Inhibit_up = Inhibit_up[!duplicated(Inhibit_up$ID),]; rownames(Inhibit_up) = Inhibit_up$ID

OPC_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "OPC"), ]
OPC_up = OPC_up[!duplicated(OPC_up$ID),]; rownames(OPC_up) = OPC_up$ID

VLMC_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "VLMC"), ]
VLMC_up = VLMC_up[!duplicated(VLMC_up$ID),]; rownames(VLMC_up) = VLMC_up$ID

Endothelial_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Endothelial"), ]
Endothelial_up = Endothelial_up[!duplicated(Endothelial_up$ID),]; rownames(Endothelial_up) = Endothelial_up$ID

pid = c("GO:0016579", "GO:0080135", "GO:0000902", "GO:0006635", "GO:1990778", "GO:0009081", "GO:0002253", "GO:0045321",
        "GO:0099536", "GO:0006024", "GO:0016358", "GO:0034329", "GO:0090083", "GO:0070482", "GO:0042026", "GO:0006458")

Excit_up = Excit_up[pid, ]
Astro_up = Astro_up[pid, ]
Oligo_up = Oligo_up[pid, ]
Microglia_up = Microglia_up[pid, ]
Inhibit_up = Inhibit_up[pid, ]
OPC_up = OPC_up[pid, ]
VLMC_up = VLMC_up[pid, ]
Endothelial_up = Endothelial_up[pid, ]

pdat = data.frame(row.names = pid,
                  Excit_up = Excit_up$score,
                  Astro_up = Astro_up$score,
                  Oligo_up = Oligo_up$score,
                  Microglia_up = Microglia_up$score,
                  Inhibit_up = Inhibit_up$score,
                  OPC_up = OPC_up$score,
                  VLMC_up = VLMC_up$score,
                  Endothelial_up = Endothelial_up$score)

pdat[is.na(pdat)] = 0

pdat = as.data.frame(t(pdat))

for(i in 1:ncol(pdat)){
  pdat[,i] = as.numeric(pdat[,i])
}

pdat$group = rownames(pdat)

pdat = pdat[,c(17,1:16)]
pdat$group = factor(pdat$group, levels = c("Astro_up", "Endothelial_up", "Excit_up", "Inhibit_up", "Microglia_up", "Oligo_up", "OPC_up", "VLMC_up"))

ggradar(pdat,
        grid.max = max(pdat[,-1]),                 # 设置坐标轴的最大值
        grid.mid = max(pdat[,-1])/2,               # 设置坐标轴的中间值
        grid.min = 0,                            # 设置坐标轴的最小值
        grid.label.size = 4,                     # 坐标轴百分比标签大小
        axis.label.size = 5,                     # 组名标签字体大小
        background.circle.colour = "white",      # 设置背景颜色
        group.point.size = 2,                    # 点大小
        group.line.width = 1,                    # 线条粗细
        plot.legend = T,                         # 是否显示图例
        legend.position = "right",               # 图例位置"top", "right", "bottom", "left"
        legend.title = "",                       # 图例标题
        legend.text.size = 10,                   # 图例文字大小
        plot.title   = "Title",                  # 标题名称
        plot.extent.x.sf = 1.2,                  # 设置图片横向延伸空间，防止外圈文字显示不全
        plot.extent.y.sf = 1.2,                  # 设置图片纵向延伸空间，防止外圈文字显示不全
) + 
  scale_color_manual(values = col_vec)

ggsave(filename = "./4-celltype/9-go-cluster.pdf", width = 9, height = 4.5)


# 可视化2
enrich_df = read.csv(file = "./4-celltype/8-enrich_outputs_go/1-GO_enrichment_all.csv")

Excit_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Excit"), ]
Excit_down = Excit_down[!duplicated(Excit_down$ID),]; rownames(Excit_down) = Excit_down$ID

Astro_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Astro"), ]
Astro_down = Astro_down[!duplicated(Astro_down$ID),]; rownames(Astro_down) = Astro_down$ID

Oligo_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Oligo"), ]
Oligo_down = Oligo_down[!duplicated(Oligo_down$ID),]; rownames(Oligo_down) = Oligo_down$ID

Microglia_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Microglia"), ]
Microglia_down = Microglia_down[!duplicated(Microglia_down$ID),]; rownames(Microglia_down) = Microglia_down$ID

Inhibit_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Inhibit"), ]
Inhibit_down = Inhibit_down[!duplicated(Inhibit_down$ID),]; rownames(Inhibit_down) = Inhibit_down$ID

OPC_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "OPC"), ]
OPC_down = OPC_down[!duplicated(OPC_down$ID),]; rownames(OPC_down) = OPC_down$ID

VLMC_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "VLMC"), ]
VLMC_down = VLMC_down[!duplicated(VLMC_down$ID),]; rownames(VLMC_down) = VLMC_down$ID

Endothelial_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Endothelial"), ]
Endothelial_down = Endothelial_down[!duplicated(Endothelial_down$ID),]; rownames(Endothelial_down) = Endothelial_down$ID

pid = c("GO:0048812", "GO:0007268", "GO:0002485", "GO:1902855", "GO:0051604", "GO:0098930", "GO:1902600", "GO:0042773",
        "GO:0006897", "GO:0016049", "GO:0016236", "GO:0061024", "GO:0046034", "GO:0006735", "GO:0071711", "GO:0031589")

Excit_down = Excit_down[pid, ]
Astro_down = Astro_down[pid, ]
Oligo_down = Oligo_down[pid, ]
Microglia_down = Microglia_down[pid, ]
Inhibit_down = Inhibit_down[pid, ]
OPC_down = OPC_down[pid, ]
VLMC_down = VLMC_down[pid, ]
Endothelial_down = Endothelial_down[pid, ]

pdat = data.frame(row.names = pid,
                  Excit_down = Excit_down$score,
                  Astro_down = Astro_down$score,
                  Oligo_down = Oligo_down$score,
                  Microglia_down = Microglia_down$score,
                  Inhibit_down = Inhibit_down$score,
                  OPC_down = OPC_down$score,
                  VLMC_down = VLMC_down$score,
                  Endothelial_down = Endothelial_down$score)

pdat[is.na(pdat)] = 0

pdat = as.data.frame(t(pdat))

for(i in 1:ncol(pdat)){
  pdat[,i] = as.numeric(pdat[,i])
}

pdat$group = rownames(pdat)

pdat = pdat[,c(17,1:16)]
pdat$group = factor(pdat$group, levels = c("Astro_down", "Endothelial_down", "Excit_down", "Inhibit_down", "Microglia_down", "Oligo_down", "OPC_down", "VLMC_down"))

ggradar(pdat,
        grid.max = max(pdat[,-1]),                 # 设置坐标轴的最大值
        grid.mid = max(pdat[,-1])/2,               # 设置坐标轴的中间值
        grid.min = 0,                            # 设置坐标轴的最小值
        grid.label.size = 4,                     # 坐标轴百分比标签大小
        axis.label.size = 5,                     # 组名标签字体大小
        background.circle.colour = "white",      # 设置背景颜色
        group.point.size = 2,                    # 点大小
        group.line.width = 1,                    # 线条粗细
        plot.legend = T,                         # 是否显示图例
        legend.position = "right",               # 图例位置"top", "right", "bottom", "left"
        legend.title = "",                       # 图例标题
        legend.text.size = 10,                   # 图例文字大小
        plot.title   = "Title",                  # 标题名称
        plot.extent.x.sf = 1.2,                  # 设置图片横向延伸空间，防止外圈文字显示不全
        plot.extent.y.sf = 1.2,                  # 设置图片纵向延伸空间，防止外圈文字显示不全
) + 
  scale_color_manual(values = col_vec)

ggsave(filename = "./4-celltype/9-go-cluster2.pdf", width = 9, height = 4.5)


