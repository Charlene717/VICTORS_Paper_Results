# Seurat V5 Azimuth annotation
# https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html
# https://azimuth.hubmapconsortium.org/
#   Azimuth outputs Prediction scores, which range from 0 to 1 and reflect the confidence associated with each annotation. However, Azimuth does not provide a recommended threshold. Can I set it at 0.8 for our tests?

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

# pkgbuild::check_build_tools(debug = TRUE)
# Download and Install Rtools: https://cran.r-project.org/bin/windows/Rtools/
# if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild")
## BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "TFBSTools"))
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if(!require("Azimuth")) devtools::install_github("satijalab/azimuth", "seurat5"); library(Azimuth)

## GitHub
if(!require("remotes")) install.packages("remotes"); library(remotes)
if(!require("SeuratData")) remotes::install_github("satijalab/seurat-data"); library(SeuratData)


# Install the PBMC systematic comparative analyis (pmbcsca) dataset
InstallData("pbmcsca")

# returns a Seurat object named pbmcsca
pbmcsca <- LoadData("pbmcsca")

# The RunAzimuth function can take a Seurat object as input
pbmcsca <- RunAzimuth(pbmcsca, reference = "pbmcref")


p1 <- DimPlot(pbmcsca, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(pbmcsca, group.by = "Method")
p1 + p2


#######################################################
# ## Error
#
# #### Load Data ####
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")
#
# seuratObject_Sample <- pbmc.rna
# seuratObject_Ref <- pbmc.rna
# seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]
#
#
# # Run Azimuth for cell type annotation
# seuratObject_Sample <- RunAzimuth(
#   query = seuratObject_Sample,
#   reference = list(
#     ref = seuratObject_Ref,
#     celltype_column = "Actual_Cell_Type"
#   ))
#
# seuratObject_Sample <- RunAzimuth(
#   query = seuratObject_Sample,
#   reference = seuratObject_Ref,
#   normalization.method = "SCT",
#   reference.assay = "refAssay",
#   query.assay = "RNA",
#   reduction = "spca",
#   dims = 1:30
# )
#
