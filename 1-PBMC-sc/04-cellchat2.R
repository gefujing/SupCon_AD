
# 设置环境
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(ggpubr)
library(CellChat)
library(patchwork)
library(ggsci)
load(file = "./4-cellchat/hc.cellchat.Rdata")
load(file = "./4-cellchat/ad.cellchat.Rdata")
load(file = "./1-qc/mycol.Rdata")

# 合并
object.list = list(hc = hc.cellchat, ad = ad.cellchat)
cellchat = mergeCellChat(object.list, add.names = names(object.list))
cellchat

# 比较细胞通讯相互作用的总数和相互作用的强度
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) + scale_fill_manual(values = col_vec)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") + scale_fill_manual(values = col_vec)
gg1 + gg2
ggsave(filename = "./4-cellchat/1-compare.pdf", width = 3.5, height = 2.6)

# 不同细胞群间相互作用次数或相互作用强度的差异
## 以circle plot的形式展示第二个组别中相较于第一个组别细胞通讯发生的变化，红色为上调蓝色为下调
pdf("./4-cellchat/2-net_number_strength.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

## 热图更详细地显示相互作用的差异数量或相互作用强度
gg1 <- netVisual_heatmap(cellchat,comparison = c(2,1))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(2,1))
gg1 + gg2
ggsave(filename = "./4-cellchat/3-compare-heatmap.pdf", width = 5.9, height = 3)

## 分开展示各组的通讯情况
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("./4-cellchat/4-net_number_strength-sigle.pdf")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

# 不同细胞间相互作用次数或相互作用强度的差异（以Fib和B细胞为例）
group.cellType = c("Bcell", "Monocyte/Macrophage", "NKcell", "NKTcell", "ProNKcell", "Tcell")

object.list2 = lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat = mergeCellChat(object.list2, add.names = names(object.list2))

## 各组别的通讯情况
weight.max <- getMaxWeight(object.list2, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
pdf("./4-cellchat/5-net_number_strength-immue.pdf")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list2)) {
  netVisual_circle(object.list2[[i]]@net$count.merged,
                   weight.scale = T, label.edge= T,
                   edge.weight.max = weight.max[3],
                   edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

## 差异性circle plot通讯情况
pdf("./4-cellchat/6-net_number_compare-immue.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged",comparison = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged",comparison = c(1,2))
dev.off()

# 比较不同组别中的主要sources和targets
num.link = sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax = c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg = list()
for (i in 1:length(object.list)) {
  object.list[[i]] = netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] = netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) + xlim(c(0,2.5)) + ylim(c(0,2.5))
}
patchwork::wrap_plots(plots = gg)

## 具体信号变化
gg1 = netAnalysis_signalingChanges_scatter(object.list, idents.use = "Bcell")
gg2 = netAnalysis_signalingChanges_scatter(object.list, idents.use = "Monocyte/Macrophage")
gg3 = netAnalysis_signalingChanges_scatter(object.list, idents.use = "Tcell")
gg4 = netAnalysis_signalingChanges_scatter(object.list, idents.use = "NKcell")
patchwork::wrap_plots(plots = list(gg1, gg2, gg3, gg4))

# 比较各个信号通路的总体信息流
gg1 = rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(1,2)) + scale_fill_aaas()
gg2 = rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(1,2)) + scale_fill_aaas()
gg1 + gg2
ggsave(filename = "./4-cellchat/9-compare-details.pdf", width = 7, height = 4.7)

# 比较与每个细胞群相关的传出outgoing（或传入incoming）信号
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

## outgoing
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 9)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 9)
draw(ht1 + ht2)

## incoming
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 9, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 9, color.heatmap = "GnBu")
draw(ht1 + ht2)

## all
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all",  signaling = pathway.union, title = names(object.list)[i],  width = 5, height = 9, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1],  width = 5, height = 9, color.heatmap = "OrRd")
draw(ht1 + ht2)


# 识别上调和下调的信号配体受体对
netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = c(1:6),  comparison = c(1,2), angle.x = 90)

## 分别展示上调下调通路
gg1 = netVisual_bubble(cellchat, sources.use = c(2,10), targets.use = c(1,5,6,9),  comparison = c(2), max.dataset = 2, title.name = "Increased signaling in tumor", angle.x = 90, remove.isolate = T)
gg2 = netVisual_bubble(cellchat, sources.use = c(2,10), targets.use = c(1,5,6,9),  comparison = c(1), max.dataset = 1, title.name = "Decreased signaling in tumor", angle.x = 90, remove.isolate = T)
gg1 + gg2

## heatmap1
pathways.show = c("CXCL") 
par(mfrow = c(2,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
