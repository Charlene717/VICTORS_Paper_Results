# PMID32567229 scClassify
# https://sydneybiox.github.io/scClassify/articles/pretrainedModel.html
# https://www.embopress.org/doi/full/10.15252/msb.20199389
# Use scClassify’s default threshold of 0.7 to evaluate the scenario where scClassify is set to force complete annotation.

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

## BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("S4Vectors", "hopach", "limma"))
# BiocManager::install(c("Rhdf5lib", "rhdf5filters", "rhdf5", "sparseMatrixStats", "graph", "HDF5Array", "DelayedMatrixStats", "GSEABase", "Cepo"))
if(!require("SingleCellExperiment")) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if(!require("scClassify")) devtools::install_github("SydneyBioX/scClassify"); library(scClassify)


#### Load Data ####
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")

seuratObject_Sample <- pbmc.rna
seuratObject_Ref <- pbmc.rna


# Convert Seurat objects to SingleCellExperiment
ref_sce <- as.SingleCellExperiment(seuratObject_Ref)
sample_sce <- as.SingleCellExperiment(seuratObject_Sample)


# Prepare Reference Data
ref_sce$celltypes <- seuratObject_Ref@meta.data[["seurat_annotations"]]


# Run scClassify
result <- scClassify(
  exprsMat_train = as.matrix(counts(ref_sce)),
  cellTypes_train = ref_sce$celltypes,
  exprsMat_test = as.matrix(counts(sample_sce))
)

result.df <- result[["testRes"]][["test"]][["pearson_WKNN_limma"]][["predLabelMat"]] %>% as.data.frame()


# Add the predicted cell types to the Seurat object
sample_sce$predicted_celltype <- result.df[,ncol(result.df)]

# Update Seurat object
seuratObject_Sample$predicted_celltype <- sample_sce$predicted_celltype

# Check results
unique(seuratObject_Sample$predicted_celltype)



plot_scClassify <- DimPlot(seuratObject_Sample,group.by = "predicted_celltype", reduction = "umap")
plot_seurat <- DimPlot(seuratObject_Sample,group.by = "seurat_annotations", reduction = "umap")
# plot_scClassify_All <- DimPlot(seuratObject_Sample,group.by = "predicted_celltype_All", reduction = "umap")

plot_scClassify + plot_seurat
# plot_scClassify_All + plot_seurat





