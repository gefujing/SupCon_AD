
## gsva
rm(list=ls()) 
library(GSVA)
library(limma)
library(tidyverse)
library(org.Hs.eg.db)
library(ggrepel)
library(enrichplot)
library(msigdbr)

options(stringsAsFactors = F) 

load(file = "./2-Batch/ad.merge.Rdata")
load(file = "./1-Cleandata/mycol.Rdata")

### 准备基因集
human_KEGG = clusterProfiler::read.gmt("./0-Rawdata/c5.go.bp.v2025.1.Hs.symbols.gmt")
human_KEGG$term = as.character(human_KEGG$term)
human_KEGG  = human_KEGG %>% split(x = .$gene, f = .$term) #后续gsva要求是list，所以将他转化为list


## 准备表达矩阵和临床信息
X = as.matrix(exp)
gsvaPar = gsvaParam(X, human_KEGG, kcdf = "Gaussian")
gsva.es = gsva(gsvaPar)
gsva.es[1:4, 1:4]


## limma差异分析
group = ph$group

design = model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(gsva.es)
contrast.matrix = makeContrasts(contrasts = paste0("AD", '-', "HC"),  #"exp/ctrl"
                                levels = design)

fit1 = lmFit(gsva.es, design)                 #拟合模型
fit2 = contrasts.fit(fit1, contrast.matrix) #统计检验
efit = eBayes(fit2)                         #修正

summary(decideTests(efit, lfc=1, p.value = 1)) #统计查看差异结果
tempOutput = topTable(efit, coef = paste0("AD", '-', "HC"), n = Inf)
degs = na.omit(tempOutput) 
degs$pathways = rownames(degs)
write.csv(degs, file = "./3-degs/7-gsva_go_degs.results.csv")


## 可视化
pdata = degs
pdata$ID = str_replace_all(string = rownames(pdata), pattern = "_", replacement = " ") %>% str_to_title()
pdata$ID = str_remove_all(string = pdata$ID, pattern = "Gobp ") %>% str_to_title()
# pdata$ID = str_remove(string = pdata$ID, pattern = " ") %>% str_to_title()
pdata = pdata[order(x = pdata$t, decreasing = T),]
pdata$rank = 1:nrow(pdata)

up = c("Positive Regulation Of Macrophage Differentiation", "Neutrophil Mediated Cytotoxicity",
       "Leukocyte Mediated Immunity", "Alpha Linolenic Acid Metabolic Process", 
       "T Cell Migration", "Regulation Of Neutrophil Migration",
       "Myeloid Leukocyte Mediated Immunity", "Myeloid Progenitor Cell Differentiation",
       "Regulation Of Mononuclear Cell Migration", "Monosaccharide Metabolic Process")

normal = pdata$ID[sample(c(1600:7006),10)]
down = c(pdata$ID[sample(c(7100:7579),7)], 
         "Protein K48 Linked Ubiquitination", "Protein Branched Polyubiquitination", "Ubiquitin Dependent Protein Catabolic Process Via The C End Degron Rule Pathway")
pp = c(up, normal, down)

rownames(pdata) = pdata$ID
pdata = pdata[pp, ]

pdata$group = ifelse(pdata$t > 2, "UP", ifelse(pdata$t < -2, "DOWN", "STABLE"))
pdata = pdata[order(x = pdata$t, decreasing = F),]
pdata$ID = factor(pdata$ID, levels = pdata$ID)


ggplot(data = pdata, aes(x = ID, y = t, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c("UP" = "#36638a", "STABLE" = "#cccccc", "DOWN" = "#7bcd7b")) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "./3-degs/7-ad-GSVA.pdf", width = 6.5, height = 7)

