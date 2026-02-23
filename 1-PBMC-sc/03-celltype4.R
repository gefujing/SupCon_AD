
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
library(scop)

# 设置数据
load(file = "./3-celltype/sce.all.filt.harmony.celltype.Rdata")
load(file = "./1-qc/mycol.Rdata")
table(sce.all.filt.harmony$celltype)

sce.all.filt.harmony@meta.data$celltype = factor(sce.all.filt.harmony@meta.data$celltype, 
                                                 levels = c("NKTcell", "Tcell", "Bcell", "Monocyte/Macrophage", "ProNKcell", "NKcell"))

# 差异基因分析
p_adj_cutoff = 0.05
min_pct = 0.10
logfc_thresh = 0
test_use = "wilcox"
min_cells_per_group = 1
assay_to_use = DefaultAssay(sce.all.filt.harmony)

obj = sce.all.filt.harmony

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

# 合并汇总
deg_summary <- bind_rows(summary_list) %>%
  arrange(desc(n_sig_total))

# 查看与保存
print(deg_summary)

# 可选：写出结果
write.csv(deg_summary, "./3-celltype/2-DEG_summary_by_celltype_AD_vs_HC.csv", row.names = FALSE)

for (ct in names(markers_all)) {
  fn <- paste0("./3-celltype/3-DEGs_", gsub("[^A-Za-z0-9_]+","_", ct), "_AD_vs_HC.csv")
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

ggsave(filename = "./3-celltype/7-degs-summary.pdf", width = 5.5, height = 4)


# 可视化2
files = list.files(
  path = "./3-celltype/7-DEGs",                 # 指定文件夹
  pattern = "^3-DEGs_.*_AD_vs_HC\\.csv$",      # 只匹配文件名
  full.names = TRUE
)


read_one <- function(fp) {
  ct <- fp %>%
    basename() %>%
    sub("^3-DEGs_(.*)_AD_vs_HC\\.csv$", "\\1", .)
  
  df <- read.csv(fp, check.names = FALSE) %>%
    as_tibble()
  
  # 兼容不同Seurat版本的列名
  if (!"avg_log2FC" %in% names(df)) {
    if ("avg_logFC" %in% names(df)) df$avg_log2FC <- df$avg_logFC
    else stop("缺少 avg_log2FC / avg_logFC: ", fp)
  }
  if (!"p_val_adj" %in% names(df)) {
    # 有些导出可能叫 p_val_adj；若没有，多半文件不对
    stop("缺少 p_val_adj: ", fp)
  }
  
  # 若没带基因列名，保证有 gene 列
  if (!"gene" %in% names(df)) {
    if (!is.null(rownames(df)) && all(rownames(df) != "")) {
      df <- df %>% rownames_to_column("gene")
    } else if ("rowname" %in% names(df)) {
      df <- df %>% rename(gene = rowname)
    } else {
      warning("未找到基因名，将用顺序代替: ", fp)
      df$gene <- paste0("gene_", seq_len(nrow(df)))
    }
  }
  
  # direction 若不存在则计算；存在则以文件为准
  if (!"direction" %in% names(df)) {
    padj_cutoff <- 0.05
    df <- df %>%
      mutate(signif = p_val_adj < padj_cutoff,
             direction = case_when(
               signif & avg_log2FC > 0  ~ "up_in_AD",
               signif & avg_log2FC < 0  ~ "down_in_AD",
               TRUE                     ~ "ns"
             ))
  }
  
  df %>%
    mutate(celltype = ct,
           neglog10_p = -log10(p_val_adj + 1e-300))
}

all_deg <- map_dfr(files, read_one)


# 保证方向为因子 & 颜色顺序
all_deg <- all_deg %>%
  mutate(direction = factor(direction, levels = c("down_in_AD", "ns", "up_in_AD")),
         celltype  = factor(celltype, levels = sort(unique(celltype))))


# ---- 图1：分面火山图（每个 celltype 一格）----
padj_line <- -log10(0.05)

# 每个 celltype 挑选上/下调 -log10(p) 最大的若干基因做标注（可调 n_label）
n_label <- 3
to_label <- all_deg %>%
  filter(direction != "ns") %>%
  group_by(celltype, direction) %>%
  slice_max(order_by = neglog10_p, n = n_label, with_ties = FALSE) %>%
  ungroup()

p_volcano <- ggplot(all_deg, aes(avg_log2FC, neglog10_p, color = direction)) +
  geom_point(size = 0.7, alpha = 0.7) +
  geom_hline(yintercept = padj_line, linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ celltype, scales = "free") +
  labs(x = "avg_log2FC (AD vs HC)", y = "-log10(adj p)") +
  scale_color_manual(values = c("down_in_AD" = "#2B6CB0", "ns" = "grey80", "up_in_AD" = "#C53030")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

# 标注基因名
p_volcano <- p_volcano +
  geom_text_repel(data = to_label,
                  aes(label = gene),
                  size = 2.6, max.overlaps = 15, min.segment.length = 0,
                  box.padding = 0.25, point.padding = 0.15, seed = 1,
                  show.legend = FALSE)

ggsave(filename = "./3-celltype/7-degs-split.pdf", width = 7, height = 5)

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
outdir <- "./3-celltype/8-enrich_outputs_go"
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
Bcell_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Bcell"), ]
Bcell_up = Bcell_up[!duplicated(Bcell_up$ID),]; rownames(Bcell_up) = Bcell_up$ID

Macrophage_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Monocyte/Macrophage"), ]
Macrophage_up = Macrophage_up[!duplicated(Macrophage_up$ID),]; rownames(Macrophage_up) = Macrophage_up$ID

NKcell_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "NKcell"), ]
NKcell_up = NKcell_up[!duplicated(NKcell_up$ID),]; rownames(NKcell_up) = NKcell_up$ID

NKTcell_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "NKTcell"), ]
NKTcell_up = NKTcell_up[!duplicated(NKTcell_up$ID),]; rownames(NKTcell_up) = NKTcell_up$ID

Tcell_up = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Tcell"), ]
Tcell_up = Tcell_up[!duplicated(Tcell_up$ID),]; rownames(Tcell_up) = Tcell_up$ID


pid = c("GO:0050853", "GO:0002250", "GO:0006886", "GO:0019882", "GO:0002521", 
        "GO:1903131", "GO:0046631", "GO:0002253", "GO:0046649", "GO:0042110")

Bcell_up = Bcell_up[pid, ]
Macrophage_up = Macrophage_up[pid, ]
NKcell_up = NKcell_up[pid, ]
NKTcell_up = NKTcell_up[pid, ]
Tcell_up = Tcell_up[pid, ]


pdat = data.frame(row.names = pid,
                  Bcell_up = Bcell_up$score,
                  Macrophage_up = Macrophage_up$score,
                  NKcell_up = NKcell_up$score,
                  NKTcell_up = NKTcell_up$score,
                  Tcell_up = Tcell_up$score)

pdat[is.na(pdat)] = 0

pdat = as.data.frame(t(pdat))

for(i in 1:ncol(pdat)){
  pdat[,i] = as.numeric(pdat[,i])
}

pdat$group = rownames(pdat)

pdat = pdat[,c(11,1:10)]
pdat$group = factor(pdat$group, levels = c("Macrophage_up", "Bcell_up", "Tcell_up", "NKTcell_up", "NKcell_up"))

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

ggsave(filename = "./3-celltype/9-go-cluster.pdf", width = 9, height = 4.5)


# 可视化2
enrich_df = read.csv(file = "./3-celltype/8-enrich_outputs_go/1-GO_enrichment_all.csv")

Bcell_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Bcell"), ]
Bcell_down = Bcell_down[!duplicated(Bcell_down$ID),]; rownames(Bcell_down) = Bcell_down$ID

Macrophage_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Monocyte/Macrophage"), ]
Macrophage_down = Macrophage_down[!duplicated(Macrophage_down$ID),]; rownames(Macrophage_down) = Macrophage_down$ID

NKcell_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "NKcell"), ]
NKcell_down = NKcell_down[!duplicated(NKcell_down$ID),]; rownames(NKcell_down) = NKcell_down$ID

NKTcell_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "NKTcell"), ]
NKTcell_down = NKTcell_down[!duplicated(NKTcell_down$ID),]; rownames(NKTcell_down) = NKTcell_down$ID

Tcell_down = enrich_df[(enrich_df$direction == "Down in AD") & (enrich_df$celltype == "Tcell"), ]
Tcell_down = Tcell_down[!duplicated(Tcell_down$ID),]; rownames(Tcell_down) = Tcell_down$ID


pid = c("GO:0002252", "GO:0090135", "GO:0006914", "GO:0010498", "GO:0035914", 
        "GO:0061535", "GO:0042110", "GO:0050852", "GO:0016032", "GO:0042093")

Bcell_down = Bcell_down[pid, ]
Macrophage_down = Macrophage_down[pid, ]
NKcell_down = NKcell_down[pid, ]
NKTcell_down = NKTcell_down[pid, ]
Tcell_down = Tcell_down[pid, ]

pdat = data.frame(row.names = pid,
                  Bcell_down = Bcell_down$score,
                  Macrophage_down = Macrophage_down$score,
                  NKcell_down = NKcell_down$score,
                  NKTcell_down = NKTcell_down$score,
                  Tcell_down = Tcell_down$score)

pdat[is.na(pdat)] = 0

pdat = as.data.frame(t(pdat))

for(i in 1:ncol(pdat)){
  pdat[,i] = as.numeric(pdat[,i])
}

pdat$group = rownames(pdat)

pdat = pdat[,c(11,1:10)]
pdat$group = factor(pdat$group, levels = c("Macrophage_down", "Bcell_down", "Tcell_down", "NKTcell_down", "NKcell_down"))

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

ggsave(filename = "./3-celltype/9-go-cluster2.pdf", width = 9, height = 4.5)




