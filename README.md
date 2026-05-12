
# SupCon_AD

This repository contains analysis code for the manuscript "Supervised contrastive learning and experimental validation identify QPCT-expressing monocytes implicated in Alzheimer’s disease".

## Data sources

This study uses publicly available datasets:
- PBMC scRNA-seq: GSE226602
- Whole-blood RNA-seq: ADNI and GSE140829
- Human brain snRNA-seq: GSE157827, GSE167494 and GSE174367
- Spatial transcriptomics: GSE220442

ADNI data are available through the ADNI data portal subject to the ADNI data use agreement. No new sequencing data were generated in this study.

## Repository structure

- 1-PBMC-sc: PBMC single-cell RNA-seq analysis
- 2-PBMC-bulk: whole-blood bulk RNA-seq analysis
- 3-Model: supervised contrastive learning model, conventional machine-learning comparison and minimal-panel classifier
- 4-Brain-sn: brain single-nucleus RNA-seq and spatial transcriptomic analysis

## Model reproduction

The supervised contrastive learning model is implemented in the 3-Model directory. The main notebooks include:
- ssl_model.ipynb: supervised contrastive learning model training and evaluation
- m10_model.ipynb: consensus-gene panel modeling and validation

Input files should be generated from the public datasets listed above following the preprocessing steps described in the Methods section of the manuscript. The model uses a 421-gene consensus feature set derived from the intersection of scRNA-seq-defined myeloid markers and bulk AD-associated differentially expressed genes.

## Software

Analyses were performed using R and Python. Major packages include Seurat, Harmony, miloR, limma, clusterProfiler, GSVA, scanpy, cell2location, PyTorch, scikit-learn, SHAP and LIME.
