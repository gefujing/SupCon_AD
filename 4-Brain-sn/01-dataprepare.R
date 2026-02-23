
# set environment
rm(list = ls()) 
options(stringsAsFactors = F) 

library(Seurat)
library(dplyr)
library(stringr)
library(data.table)
library(tidyverse)
library(ggsci)





# process GSE157827
meta = read.csv(file = "./0-Rawdata/GSE157827/GSE157827_meta.csv")
assays = dir("./0-Rawdata/GSE157827/Rawcounts/")
dir = paste0("./0-Rawdata/GSE157827/Rawcounts/", assays)

samples_name = assays
scRNAlist = list()

for(i in 1:length(dir)){
  counts = Read10X(data.dir = dir[i])
  scRNAlist[[i]] = CreateSeuratObject(counts, project = samples_name[i])
  scRNAlist[[i]] = RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i]) 
}

names(scRNAlist) = samples_name

## combine
sce.GSE157827 = merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
sce.GSE157827

# get metadata
md0 = sce.GSE157827@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(orig.ident = as.character(orig.ident))

md0 = md0[,c(1:4)]

md_joined = md0 %>%
  left_join(meta, by = c("orig.ident" = "geo_accession"))

sce.GSE157827@meta.data = md_joined %>%
  column_to_rownames("cell_id")

head(sce.GSE157827@meta.data)







# process GSE167494
meta = read.csv(file = "./0-Rawdata/GSE167494/GSE167494_meta.csv")
assays = dir("./0-Rawdata/GSE167494/Rawcounts/")
dir = paste0("./0-Rawdata/GSE167494/Rawcounts/", assays)

samples_name = assays
scRNAlist = list()

for(i in 1:length(dir)){
  counts = Read10X(data.dir = dir[i])
  scRNAlist[[i]] = CreateSeuratObject(counts, project = samples_name[i])
  scRNAlist[[i]] = RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i]) 
}

names(scRNAlist) = samples_name

## combine
sce.GSE167494 = merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
sce.GSE167494

# get metadata
md0 = sce.GSE167494@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(orig.ident = as.character(orig.ident))

md0 = md0[,c(1:4)]

md_joined = md0 %>%
  left_join(meta, by = c("orig.ident" = "geo_accession"))

sce.GSE167494@meta.data = md_joined %>%
  column_to_rownames("cell_id")
head(sce.GSE167494@meta.data)









# process GSE174367
meta = read.csv(file = "./0-Rawdata/GSE174367/GSE174367_meta.csv")
cellmeta = read.csv(file = "./0-Rawdata/GSE174367/GSE174367_snRNA-seq_cell_meta.csv")
sce.GSE174367 = Read10X_h5(filename = "./0-Rawdata/GSE174367/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
sce.GSE174367 = CreateSeuratObject(sce.GSE174367)

# cell integration
sce.GSE174367 = subset(sce.GSE174367, cells = cellmeta$Barcode)
sample_vec = stats::setNames(cellmeta$SampleID, cellmeta$Barcode)
sce.GSE174367$orig.ident = sample_vec[Cells(sce.GSE174367)]


# get metadata
md0 = sce.GSE174367@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(orig.ident = as.character(orig.ident))

md0 = md0[,c(1:4)]

md_joined = md0 %>%
  left_join(meta, by = c("orig.ident" = "sample"))

sce.GSE174367@meta.data = md_joined %>%
  column_to_rownames("cell_id")

sce.GSE174367@meta.data = sce.GSE174367@meta.data[,c(4,2,3,1,5,6,7,8,9,10,11,12)]
colnames(sce.GSE174367@meta.data)[1] = "orig.ident"
colnames(sce.GSE174367@meta.data)[4] = "sample"
head(sce.GSE174367@meta.data)






## 合并数据集
sce.all = merge(x = sce.GSE157827, y = sce.GSE167494)
sce.all = merge(x = sce.all, y = sce.GSE174367)


## 保存数据
save(sce.GSE157827, file = "./1-datapepare/sce.GSE157827.Rdata")
save(sce.GSE167494, file = "./1-datapepare/sce.GSE167494.Rdata")
save(sce.GSE174367, file = "./1-datapepare/sce.GSE174367.Rdata")
save(sce.all, file = "./1-datapepare/sce.all.Rdata")



col_vec = c(pal_npg()(8),  
            pal_aaas()(8), 
            pal_nejm()(8), 
            pal_lancet()(8), 
            pal_jama()(7), 
            pal_bmj()(8), 
            pal_jco()(8),
            pal_ucscgb()(26))[sample(70,70)]
col_vec = unique(col_vec)
save(col_vec, file = "./1-datapepare/mycol.Rdata")
