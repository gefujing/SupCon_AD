
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
load("./3-celltype/sce.all.filt.harmony.celltype.Rdata")

# 取非转移子集
hc.sce = subset(x = sce.all.filt.harmony, Group == "HC")
rm(list = c("sce.all.filt.harmony"))

# 导入数据
hc.sce = SetIdent(hc.sce, value = "celltype")
DimPlot(hc.sce, reduction = "umap", group.by = "celltype", label = T)

# 输入数据
data.input = hc.sce@assays$RNA@layers$data
rownames(data.input) = rownames(hc.sce)
colnames(data.input) = colnames(hc.sce)

identity = data.frame(group = hc.sce$celltype, row.names = names(hc.sce$celltype)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels

# 创建cellchat对象
cellchat = createCellChat(object = data.input)
cellchat

# 加入metadata
cellchat = addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat = setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize = as.numeric(table(cellchat@idents)) # number of cells in each cell group

# 导入配体受体数据库
CellChatDB = CellChatDB.human
showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use = subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB = CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)

# 预处理
cellchat = subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat = identifyOverExpressedGenes(cellchat) #寻找高表达的基因#
cellchat = identifyOverExpressedInteractions(cellchat) #寻找高表达的通路
cellchat = projectData(cellchat, PPI.human) #投影到PPI


# 相互作用推断
cellchat = computeCommunProb(cellchat) #默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat = computeCommunProbPathway(cellchat)
cellchat = filterCommunication(cellchat, min.cells = 10)
df.netp = subsetCommunication(cellchat, slot.name = "netP")
cellchat = aggregateNet(cellchat)

# 保存数据
hc.cellchat = cellchat
save(hc.cellchat, file = "./4-cellchat/hc.cellchat.Rdata")

# 设置环境
rm(list=ls())
load("./3-celltype/sce.all.filt.harmony.celltype.Rdata")

# 取转移子集
ad.sce = subset(x = sce.all.filt.harmony, Group == "AD")
rm(list = c("sce.all.filt.harmony"))

# 导入数据
ad.sce = SetIdent(ad.sce, value = "celltype")
DimPlot(ad.sce, reduction = "umap", group.by = "celltype", label = T)

# 输入数据
data.input = ad.sce@assays$RNA@layers$data
rownames(data.input) = rownames(ad.sce)
colnames(data.input) = colnames(ad.sce)

identity = data.frame(group = ad.sce$celltype, row.names = names(ad.sce$celltype)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels

# 创建cellchat对象
cellchat = createCellChat(object = data.input)
cellchat

# 加入metadata
cellchat = addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat = setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize = as.numeric(table(cellchat@idents)) # number of cells in each cell group

# 导入配体受体数据库
CellChatDB = CellChatDB.human
showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use = subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB = CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)

# 预处理
cellchat = subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat = identifyOverExpressedGenes(cellchat) #寻找高表达的基因#
cellchat = identifyOverExpressedInteractions(cellchat) #寻找高表达的通路
cellchat = projectData(cellchat, PPI.human) #投影到PPI

# 相互作用推断
cellchat = computeCommunProb(cellchat) #默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat = computeCommunProbPathway(cellchat)
cellchat = filterCommunication(cellchat, min.cells = 10)
df.netp = subsetCommunication(cellchat, slot.name = "netP")
cellchat = aggregateNet(cellchat)

# 保存数据
ad.cellchat = cellchat
save(ad.cellchat, file = "./4-cellchat/ad.cellchat.Rdata")
