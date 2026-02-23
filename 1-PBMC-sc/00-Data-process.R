
# set environment
rm(list = ls())
library(Seurat)
library(tidyverse)
library(Matrix)

# process GSE226602
counts = readRDS(file = "./0-Rawdata/GSE226602/GSE226602_rna_raw_counts.rds")
meta = read.csv(file = "./0-Rawdata/GSE226602/GSE226602_metadata.csv")
sce.GSE226602 = CreateSeuratObject(counts = counts)
sce.GSE226602

# get metadata
md0 = sce.GSE226602@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(orig.ident = as.character(orig.ident))

md_joined = md0 %>%
  left_join(meta, by = c("orig.ident" = "Sample"))

sce.GSE226602@meta.data = md_joined %>%
  column_to_rownames("cell_id")


# process GSE181279
meta = read.csv(file = "./0-Rawdata/GSE181279/GSE181279_metadata.csv")
assays = dir("./0-Rawdata/GSE181279/")
assays = assays[-4]
dir = paste0("./0-Rawdata/GSE181279/", assays)

## read counts
samples_name = assays
scRNAlist = list()

for(i in 1:length(dir)){
  counts = Read10X(data.dir = dir[i])
  scRNAlist[[i]] = CreateSeuratObject(counts, project = samples_name[i])
  scRNAlist[[i]] = RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i]) 
}

names(scRNAlist) = samples_name

## combine
sce.GSE181279 = merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
sce.GSE181279


# get metadata
md0 = sce.GSE181279@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(orig.ident = as.character(orig.ident))

md_joined = md0 %>%
  left_join(meta, by = c("orig.ident" = "Sample"))

sce.GSE181279@meta.data = md_joined %>%
  column_to_rownames("cell_id")

# merge data
sce.all = merge(sce.GSE181279, sce.GSE226602)
sce.all = JoinLayers(sce.all)
save(sce.all, file = "./1-Dataprepare/sce.all.Rdata")


# output data
DefaultAssay(sce.all) = "RNA"
mat = GetAssayData(sce.all, slot = "counts")
stopifnot(inherits(mat, "dgCMatrix"))
stopifnot(identical(colnames(mat), rownames(sce.all@meta.data)))

# write matrix
Matrix::writeMM(mat, file = "./1-Dataprepare/matrix.mtx")
write_tsv(tibble(rownames(mat)), file = "./1-Dataprepare/features.tsv", col_names = FALSE)
write_tsv(tibble(colnames(mat)), file = "./1-Dataprepare/barcodes.tsv", col_names = FALSE)

# zip all files
R.utils::gzip("./1-Dataprepare/matrix.mtx", overwrite = TRUE) 
R.utils::gzip("./1-Dataprepare/features.tsv", overwrite = TRUE)
R.utils::gzip("./1-Dataprepare/barcodes.tsv", overwrite = TRUE) 

# write metadata
md = sce.all@meta.data |>
  rownames_to_column("barcode")
readr::write_csv(md, "./1-Dataprepare/metadata.csv.gz")





