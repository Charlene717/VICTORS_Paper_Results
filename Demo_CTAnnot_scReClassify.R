# PMID31874628 scReClassify
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6305-x
# https://bioconductor.org/packages/release/bioc/vignettes/scReClassify/inst/doc/scReClassify.html

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

# if (!("devtools" %in% rownames(installed.packages()))) install.packages("devtools"); library(devtools)
# if(!require("scReClassify")) devtools::install_github("SydneyBioX/scReClassify"); library(scReClassify)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("scReClassify")) BiocManager::install("scReClassify"); library(scReClassify)



#### Load Data ####
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")

seuratObject_Sample <- pbmc.rna
seuratObject_Ref <- pbmc.rna
