# SupCon_AD

This repository contains the analysis code for the manuscript "Supervised contrastive learning and experimental validation identify QPCT-expressing monocytes implicated in Alzheimer’s disease".

## System requirements

The code was developed and tested on a standard Linux-based computational environment. No non-standard hardware is required for running the preprocessing and downstream analysis scripts. GPU acceleration is recommended but not strictly required for training the PyTorch supervised contrastive learning model.

Major software dependencies include:
- Python 3.9 or later
- R 4.2 or later
- PyTorch
- scikit-learn
- pandas
- numpy
- scipy
- matplotlib
- seaborn
- scanpy
- cell2location
- shap
- lime
- Seurat
- Harmony
- miloR
- limma
- clusterProfiler
- GSVA

## Installation

Clone the repository:

git clone https://github.com/gefujing/SupCon_AD.git
cd SupCon_AD

Install Python dependencies using conda or pip according to the package requirements of the analysis notebooks. R dependencies should be installed from CRAN, Bioconductor or GitHub as appropriate.

Typical installation time on a normal desktop computer is approximately 30-60 minutes, depending on the operating system and package cache status.

## Data sources

This study uses publicly available datasets:
- PBMC scRNA-seq: GSE226602
- Whole-blood RNA-seq: ADNI and GSE140829
- Human brain snRNA-seq: GSE157827, GSE167494 and GSE174367
- Spatial transcriptomics: GSE220442

ADNI data are available through the ADNI data portal and are subject to the ADNI data use agreement. No new sequencing data were generated in this study.

## Repository structure

- 1-PBMC-sc: PBMC single-cell RNA-seq analysis
- 2-PBMC-bulk: whole-blood bulk RNA-seq analysis
- 3-Model: supervised contrastive learning model, conventional machine-learning comparison and minimal-panel classifier
- 4-Brain-sn: brain single-nucleus RNA-seq and spatial transcriptomic analysis

## Running the model

The supervised contrastive learning model is implemented in the 3-Model directory.

Main notebooks:
- ssl_model.ipynb: supervised contrastive learning model training and evaluation
- m10_model.ipynb: consensus-gene panel modeling and validation

The model uses a 421-gene consensus feature set derived from the intersection of scRNA-seq-defined myeloid markers and bulk AD-associated differentially expressed genes. Input expression matrices should be log1p-transformed and z-scored per gene, as described in the manuscript Methods.

Expected outputs include trained model checkpoints, prediction probabilities, cross-validation performance metrics, ROC/PR curves, confusion matrices and feature-importance rankings.

## Demo dataset

A small simulated demo dataset can be provided in the 0-test-data directory to demonstrate code execution without redistributing controlled-access clinical data.

Expected runtime for the demo on a normal desktop computer is approximately 5-15 minutes. Full analyses using public single-cell, bulk and spatial transcriptomic datasets require substantially longer runtime depending on dataset size and computing resources.

## Reproduction instructions

To reproduce the main computational analyses, download the public datasets listed above, preprocess them according to the Methods section of the manuscript, and run the scripts/notebooks in the corresponding repository folders. The supervised contrastive learning and consensus-panel analyses can be reproduced using the notebooks in the 3-Model directory.

## License

The software is released under the MIT License.
