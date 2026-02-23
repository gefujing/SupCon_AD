
# 加载必要的库
rm(list = ls())
options(stringsAsFactors = F)
library(Mfuzz)
library(edgeR)
library(limma)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ClusterGVis)
library(RColorBrewer)
library(org.Hs.eg.db)
library(ggsci)
library(ggpubr)


# 数据读取
exp = fread(input = "./0-Rawdata/ADNI/adni.matrix.csv", data.table = F)

exp = exp[!duplicated(exp$GeneID), ]
exp = exp[!str_detect(string = exp$GeneID, pattern = " "), ]
rownames(exp) = exp$GeneID
exp = exp[,-1]

# match phenotype
ph = read.csv(file = "./0-Rawdata/ADNI/adni.ph.csv")
rownames(ph) = ph$PTID

kk = intersect(rownames(ph), colnames(exp))

ph = ph[kk,]
exp = exp[,kk]

# 数据排序
ph$DX.bl = factor(ph$DX.bl, levels = c("CN", "EMCI", "LMCI", "AD"))
ph = ph[order(ph$DX.bl),]

exp = exp[,rownames(ph)]
mps = as.matrix(exp)
colnames(mps) = ph$DX.bl

# 数据预处理
mps = log2(mps + 1)
mps = t(limma::avereps(t(mps) , ID = colnames(mps))) # 取平均值
mps = na.omit(mps)

ck = clusterData(exp = mps, cluster.method = "kmeans", cluster.num = 4)
visCluster(object = ck, plot.type = "line", ncol = 4)
ggsave(filename = "./4-adni-degs/1-gene-change.pdf", width = 9.5, height = 2.8)


# 聚类
genelist = ck[["wide.res"]]
write.csv(genelist, file = "./4-adni-degs/genecluster.csv")

s2e = bitr(genelist$gene, 
           fromType = "SYMBOL",
           toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)

genelist = inner_join(genelist, s2e, by=c("gene" = "SYMBOL"))

enrich_KEGG <- data.frame()
for (i in 1:4) {
  p = enrichKEGG(
    gene = genelist$ENTREZID[which(genelist$cluster == i)],
    keyType = 'kegg',
    organism = 'hsa',  
    pAdjustMethod = 'fdr',  
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.2) 
  p = setReadable(p, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  pp = data.frame(id=paste0("C",i), term=p@result[["Description"]], p.adjust=p@result[["p.adjust"]], gene=p@result[["geneID"]])
  enrich_KEGG <- rbind(enrich_KEGG, pp)
}

write.csv(enrich_KEGG, file = "./4-adni-degs/keggcluster.csv")


cg = genelist$ENTREZID[genelist$cluster == 4]
p = enrichKEGG(
  gene = cg,
  keyType = 'kegg',
  organism = 'hsa',  
  pAdjustMethod = 'fdr',  
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.2) 
p = setReadable(p, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
p = p@result


p = enrichGO(
  gene = cg,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pvalueCutoff = 1, 
  qvalueCutoff = 1) 
p = setReadable(p, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
p = p@result

write.csv(p, file = "./4-adni-degs/go4cluster.csv")


#绘图
p$ONTOLOGY = 1:nrow(p)
kkchange = p[c(3,4,7,20,26,41,56,64,110,367,401,448,463,675,685), ]

kkchange$order = factor(rev(as.integer(1:15)),
                        labels = rev(kkchange$Description))
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
  labs(title = "GO enrichment",
       x = "Count", 
       y = "Pathways")+
  theme_bw()+
  theme(text = element_text(colour = "black"))
p2
ggsave(p2, filename = "./4-adni-degs/2-KEGG.pdf.pdf", width = 6, height = 3.5)


# 组合图（热图+箱线图）
markgenes = c("APP", "PSEN1", "MME", "APOE", "CLU", "LRP1",
              "C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "C5")
visCluster(object = ck,
           plot.type = "both",
           markGenes = markgenes,
           annoTerm.data = enrich_KEGG, 
           # line.side = "left", 
           # markGenes = mark,
           show_row_dend = F)



