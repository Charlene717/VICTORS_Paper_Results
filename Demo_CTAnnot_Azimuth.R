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

# Download and Install Rtools: https://cran.r-project.org/bin/windows/Rtools/
# if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild")
# pkgbuild::check_build_tools(debug = TRUE)
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




















#### Load Data ####
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")

seuratObject_Sample <- pbmc.rna
seuratObject_Ref <- pbmc.rna


# Run Azimuth for cell type annotation
query <- RunAzimuth(
  object = seuratObject_Sample,
  reference = seuratObject_Ref,
  normalization.method = "SCT",
  reference.assay = "refAssay",
  query.assay = "RNA",
  reduction = "spca",
  dims = 1:30
)

