
# 设置环境
rm(list = ls()) 
options(stringsAsFactors = F) 

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)

# 设置数据
load(file = "./5-macsub/macsub.Rdata")
load(file = "./1-qc/mycol.Rdata")
table(macsub$celltype)

macsub = as.SingleCellExperiment(macsub)
macsub = Milo(macsub)

# Construct KNN
traj_milo = buildGraph(macsub, k = 10, d = 30)

# Defining neighbourhoods
traj_milo = makeNhoods(traj_milo, prop = 0.1, k = 10, d = 30, refined = TRUE)

# Counting cells in neighbourhoods
traj_milo = countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples = "orig.ident")

# Differential abundance testing
traj_design = data.frame(colData(traj_milo))[,c("orig.ident", "Group")]
traj_design$Group = factor(traj_design$Group, levels = c("HC", "AD"))
traj_design = distinct(traj_design)
rownames(traj_design) = traj_design$orig.ident

## Reorder rownames to match columns of nhoodCounts(milo)
traj_design = traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]
traj_milo = calcNhoodDistance(traj_milo, d=30)
rownames(traj_design) = traj_design$orig.ident
da_results = testNhoods(traj_milo, design = ~ Group, design.df = traj_design)
da_results = annotateNhoods(traj_milo, da_results, coldata_col = "celltype")

write.csv(x = da_results, file = "./5-macsub/4-miloR.csv")

##Visualize neighbourhoods displaying DA
traj_milo = buildNhoodGraph(traj_milo)

plotNhoodGraphDA(traj_milo, da_results, alpha=0.05)
ggsave(filename = "./5-macsub/10-milo.pdf", width = 6.5, height = 5)

plotDAbeeswarm(da_results, group.by = "celltype", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = 4)+
  theme_classic()+
  scale_color_gradient2(low = "navy", mid = "white", high = "#7d1c1c")
ggsave(filename = "./5-macsub/11-milo.pdf", width = 5, height = 3)


# 可视化
ggplot (da_results , aes ( logFC , - log10 ( SpatialFDR ), color = celltype)) +
  geom_point () +
  geom_hline ( yintercept  =  1 )
