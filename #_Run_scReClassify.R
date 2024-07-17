# PMID31874628 scReClassify
# https://bioconductor.org/packages/release/bioc/vignettes/scReClassify/inst/doc/scReClassify.html


# ##### Presetting ######
# rm(list = ls()) # Clean variable ##* Comment out if Run All
# memory.limit(150000)

#### Load Packages ####
# Load necessary packages
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("caret")) install.packages("caret"); library(caret)

if (!require("scReClassify")) BiocManager::install("scReClassify"); library(scReClassify)
if (!require("DT")) install.packages("DT"); library(DT)
if (!require("mclust")) install.packages("mclust"); library(mclust)
if(!require("SingleCellExperiment"))  BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
