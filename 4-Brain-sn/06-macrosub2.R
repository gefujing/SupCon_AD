
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


# 导入数据
load("./1-datapepare/mycol.Rdata")
load("./6-microsub/microsub.Rdata")




