# PMID31226206 CHETAH
# https://www.bioconductor.org/packages/devel/bioc/vignettes/CHETAH/inst/doc/CHETAH_introduction.html
# Use CHETAH’s default threshold of 0.1 to evaluate the scenario where CHETAH is set to force complete annotation.

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

## BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("SingleCellExperiment")) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
if(!require("CHETAH")) BiocManager::install("CHETAH"); library(CHETAH)


#### Load Data ####
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")

seuratObject_Sample <- pbmc.rna
seuratObject_Ref <- pbmc.rna
# seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]

# Convert Seurat objects to SingleCellExperiment
ref_sce <- as.SingleCellExperiment(seuratObject_Ref)
sample_sce <- as.SingleCellExperiment(seuratObject_Sample)


# Prepare Reference Data
ref_sce$celltypes <- seuratObject_Ref@meta.data[["seurat_annotations"]]


# Run CHETAH classifier
sample_sce <- CHETAHclassifier(input = sample_sce, ref_cells = ref_sce)
sample_sce_All <- CHETAHclassifier(input = sample_sce, ref_cells = ref_sce, thresh = 0)

# Plot classification
PlotCHETAH(sample_sce)

# Extract cell types
celltypes <- sample_sce$celltype_CHETAH
celltypes_All <- sample_sce_All$celltype_CHETAH


# Update Seurat Object
seuratObject_Sample$predicted_celltype <- celltypes
seuratObject_Sample$predicted_celltype_All <- celltypes_All



plot_CHETAH <- DimPlot(seuratObject_Sample,group.by = "predicted_celltype", reduction = "umap")
plot_seurat <- DimPlot(seuratObject_Sample,group.by = "seurat_annotations", reduction = "umap")
plot_CHETAH_All <- DimPlot(seuratObject_Sample,group.by = "predicted_celltype_All", reduction = "umap")

plot_CHETAH + plot_seurat
plot_CHETAH_All + plot_seurat





##################################################################################
# # https://www.bioconductor.org/packages/devel/bioc/vignettes/CHETAH/inst/doc/CHETAH_introduction.html
#
# #### Load Data ####
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/CHETAH_TME/CHETAH_TME_reference.Rdata")
#
# ## Make SingleCellExperiments
# reference <- SingleCellExperiment(assays = list(counts = ref_counts),
#                                   colData = DataFrame(celltypes = ref_ct))
#
# input <- SingleCellExperiment(assays = list(counts = input_counts),
#                               reducedDims = SimpleList(TSNE = input_tsne))
#
# ## Run CHETAH
# input <- CHETAHclassifier(input = input, ref_cells = reference)
#
# ## Plot the classification
# PlotCHETAH(input)
#
# ## Extract celltypes:
# celltypes <- input$celltype_CHETAH
#
