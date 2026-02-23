
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(limma)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(AnnoProbe)
library(tinyarray)
library(data.table)
library(ggrepel)
library(eulerr)
library(ggvenn)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(aPEAR)
library(GOplot)
library(GseaVis)
library(ggsci)
library(ggpubr)

load(file = "./2-Batch/ad.merge.Rdata")
load(file = "./1-Cleandata/mycol.Rdata")

# 数据预处理
grouplist = ph$group

# limma分析差异基因
design = model.matrix(~0 + grouplist)
colnames(design) = levels(grouplist)
rownames(design) = colnames(exp)

fit = lmFit(exp, design)

constrasts = paste(rev(levels(grouplist)), collapse = "-")
cont.matrix = makeContrasts(contrasts = constrasts, levels = design) 
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

deg = topTable(fit2, coef = constrasts, n = Inf)
deg = na.omit(deg)

# 差异基因统计分析
pvalue_t = 0.05
logFC_t = 0
k1 = (deg$P.Value < pvalue_t)&(deg$logFC < -logFC_t);table(k1)
k2 = (deg$P.Value < pvalue_t)&(deg$logFC > logFC_t);table(k2)
deg$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))

table(deg$change)

## 保存数据
write.csv(deg, file = "./3-degs/1-advshc.degs.csv")
save(deg, file = "./3-degs/1-advshc.degs.Rdata")

write.csv(exp, file = "./3-degs/exp_all.csv")
write.csv(ph, file = "./3-degs/ph_all.csv")

# save
cg = rownames(deg)[deg$change != "NOT"]
exp = exp[cg,]

write.csv(exp, file = "./3-degs/exp.csv")
write.csv(ph, file = "./3-degs/ph.csv")

# 可视化
pdata = deg
pdata$symbol = rownames(pdata)
for_label = pdata[pdata$symbol[1:5],]
# for_label = pdata[pdata$symbol %in% c("IL1B", "CCL3", "CCL4"),]

ggplot(data = pdata) + 
  geom_point(aes(x = logFC, y = -log10(P.Value), color = logFC, size = -log10(P.Value))) + 
  geom_text_repel(data =  for_label, aes(x = logFC, y = -log10(P.Value), label = symbol),
                  nudge_x = 0.1, nudge_y = 0.2, segment.curvature = -0.1,
                  segment.ncp = 1, direction = "y", hjust = "left", max.overlaps = 200 ) + 
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(log2(1)), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  theme_bw() + 
  ggtitle(label = "AD(248) vs. HC(509)")
ggsave(filename = "./3-degs/1-degs-advshc.pdf", width = 5, height = 3.5)


# 富集分析KEGG
deg$symbol = rownames(deg)
s2e = bitr(deg$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
deg = inner_join(deg, s2e, by=c("symbol" = "SYMBOL"))

# 输入数据
gene_diff = deg$ENTREZID[deg$change != "NOT"] 

# KEGG富集
ekk = enrichKEGG(gene = gene_diff, organism = "hsa")
ekk = setReadable(ekk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
ekk_table = ekk@result
write.csv(x = ekk_table, file = "./3-degs/2-kegg_enrich.csv")


#绘图
ekk_table$ID = 1:nrow(ekk_table)
kkchange = ekk_table[c(7,8,9,10,11,12,15,18,19,21,23,25,27,28,29,31,32,33,34,35) , ]

kkchange$order = factor(as.integer(1:20), labels = kkchange$Description)
fz = str_split(kkchange$GeneRatio, pattern = "/", simplify = T)[,1]
fm = str_split(kkchange$GeneRatio, pattern = "/", simplify = T)[,2]
kkchange$GeneRatio = as.numeric(fz)/as.numeric(fm)

fz = str_split(kkchange$BgRatio, pattern = "/", simplify = T)[,1]
fm = str_split(kkchange$BgRatio, pattern = "/", simplify = T)[,2]
kkchange$BgRatio = as.numeric(fz)/as.numeric(fm)

kkchange$EnrichFactor = kkchange$GeneRatio/kkchange$BgRatio

p2 = ggplot(kkchange, aes(x = Count, y = order, fill = pvalue)) + 
  geom_bar(stat='identity') +
  scale_fill_gradient(low = "red",high = "yellow" ) +
  labs(title = "KEGG",
       x = "Count", 
       y = "Pathways")+
  theme_classic()+
  scale_x_continuous(expand = c(0,0)) +
  coord_flip()+
  theme(text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p2
ggsave(p2, filename = "./3-degs/2-KEGG.pdf.pdf", width = 6, height = 4.5)


# GO富集分析
ego = enrichGO(gene = gene_diff, OrgDb = org.Hs.eg.db, ont = "ALL", readable = TRUE)
ego_table = ego@result
write.csv(x = ego_table, file = "./3-degs/3-go_enrich.csv")

enrichmentNetwork(ego@result[1:50,], colorBy = 'pvalue', colorType = 'pval', innerCutoff = 0.1, outerCutoff = 0.5, pCutoff = -5)
ggsave(filename = "./3-degs/3-GO1.pdf.pdf", width = 7.5, height = 4.5)

# GO可视化
ego_table$GeneRatio = 1:nrow(ego_table)
GO = ego_table[c(2,3,7,8,13,14,18,32,200,283), c(2,3,9,6)]
GO$geneID <- str_replace_all(GO$geneID,"/",",") ### 修改geneID这一列的分隔符号
names(GO)=c("ID","Term","Genes","adj_pval")
GO$Category = "BP"

genedata = data.frame(ID = deg$symbol, logFC = deg$logFC)
circ = circle_dat(GO, genedata)

gene = circ %>% group_by(term) %>% top_n(5, logFC)
chord = chord_dat(circ, genes = gene$genes)
GOChord(chord, ribbon.col = col_vec[1:10])

ggsave(filename = "./3-degs/4-GO2.pdf", width = 5.7, height = 7)


# ID转换
deg$symbol = rownames(deg)
s2e = bitr(deg$symbol, 
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)

deg = inner_join(deg, s2e, by=c("symbol" = "SYMBOL"))


## 数据构建
deg = deg[order(deg$logFC, decreasing = T),]
genelist = deg$logFC
names(genelist) = deg$ENTREZID

## 利用KEGG进行富集
gsea_go = gseGO(
  geneList = genelist, # 根据logFC排序的基因集
  OrgDb = "org.Hs.eg.db" 

)

## 导出表格
gsea_go_table = gsea_go@result
write.csv(gsea_kk_table, file = "./3-degs/4-gsea.csv")

## 可视化
volcanoGsea(data = gsea_go,
            topN = 5,
            nudge.y = c(-0.8,0.8)
)


geneSetID = c("GO:0002274", "GO:0042119", "GO:0051402", "GO:0070936", "GO:0045333", "GO:0071826", "GO:0042254")
gseaNb(object = gsea_go,
       geneSetID = geneSetID[1:3],
       addPval = T)

gseaNb(object = gsea_go,
       geneSetID = geneSetID[4:6],
       addPval = T)


# AD相关基因
# 数据读取
cg = c('ARHGAP17', 'BTK', 'NAGK', 'PAK1', 'PRR13', 'QPCT', 'RYBP', 'SLC22A4', 'TAGLN', 'TMEM154')

pdata = ph[,c(1,2)]
colnames(pdata) = c("id", "Group")
pdata[,cg] = t(exp[cg,])

pdata = pivot_longer(data = pdata, cols = cg, values_to = "Expression", names_to = "Genes")
pdata$Genes = factor(pdata$Genes, levels = cg)
pdata$Group = factor(pdata$Group, levels = c("HC", "AD"))
pdata$Expression = log2(pdata$Expression)


p2 = ggplot(pdata, aes(x = Genes, y = Expression, color = Group))+
  geom_boxplot(aes(fill=Group), alpha=0.1)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), size = 0.5)+
  scale_color_manual(values = pal_npg('nrc')(9)[c(1,2,3)])+
  scale_fill_manual(values = pal_npg('nrc')(9)[c(1,2,3)])+
  ylab("Relative expression")+
  theme_classic()

p2 <- p2 + stat_compare_means(aes(group = Group),
                              label="p.signif",
                              show.legend = F)
p2

ggsave(p2, filename = "./3-degs/8-myeloid-genes2.pdf", width = 7, height = 4)
















