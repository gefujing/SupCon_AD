

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
load("./5-neusub/neusub.Rdata")
load(file = "./1-datapepare/mycol.Rdata")

# 亚群再分析
exisub = subset(x = neusub, celltype == "Excit")
inhsub = subset(x = neusub, celltype == "Inhibit")

exideg = RunDEtest(srt = exisub, group_by = "diagnosis", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = exideg, group_by = "diagnosis")
ggsave(filename = "./5-neusub/5-exi-degs.pdf", width = 6.5, height = 4)

DEGs = exideg@tools$DEtest_diagnosis$AllMarkers_wilcox
DEGs = DEGs[with(DEGs, p_val_adj < 0.05), ]
write.csv(DEGs, file = "./5-neusub/1-exi-degs.csv")



inhdeg = RunDEtest(srt = inhsub, group_by = "diagnosis", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = inhdeg, group_by = "diagnosis")
ggsave(filename = "./5-neusub/5-inh-degs.pdf", width = 6.5, height = 4)

DEGs = inhdeg@tools$DEtest_diagnosis$AllMarkers_wilcox
DEGs = DEGs[with(DEGs, p_val_adj < 0.05), ]
write.csv(DEGs, file = "./5-neusub/1-inh-degs.csv")




