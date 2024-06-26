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

if(!require("SingleCellExperiment")) install.packages("SingleCellExperiment"); library(SingleCellExperiment)


## Make SingleCellExperiments
reference <- SingleCellExperiment(assays = list(counts = ref_counts),
                                  colData = DataFrame(celltypes = ref_ct))

input <- SingleCellExperiment(assays = list(counts = input_counts),
                              reducedDims = SimpleList(TSNE = input_tsne))

## Run CHETAH
input <- CHETAHclassifier(input = input, ref_cells = reference)

## Plot the classification
PlotCHETAH(input)

## Extract celltypes:
celltypes <- input$celltype_CHETAH

