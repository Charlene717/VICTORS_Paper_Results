##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("caret")) install.packages("caret"); library(caret)

# PMID31874628 scReClassify # https://bioconductor.org/packages/release/bioc/vignettes/scReClassify/inst/doc/scReClassify.html

# Load necessary packages
if (!require("scReClassify")) BiocManager::install("scReClassify"); library(scReClassify)
if (!require("DT")) install.packages("DT"); library(DT)
if (!require("mclust")) install.packages("mclust"); library(mclust)

# Convert Seurat objects to SingleCellExperiment
# ref_sce <- as.SingleCellExperiment(Reference_Seurat)
sample_sce <- as.SingleCellExperiment(Query_Seurat)

# Standardize the data and create 'logNorm' assay
# Log-normalize dataset if not already done
if (is.null(Reference_Seurat@assays[["RNA"]]@layers[["data"]])) { ref_sce <- logNormCounts(ref_sce)}
if (is.null(Query_Seurat@assays[["RNA"]]@layers[["data"]])) { sample_sce <- logNormCounts(sample_sce)}


# # Create 'logNorm' layer in assays
# SummarizedExperiment::assay(ref_sce, "logNorm") <- SummarizedExperiment::assay(ref_sce, "logcounts")
SummarizedExperiment::assay(sample_sce, "logNorm") <- SummarizedExperiment::assay(sample_sce, "logcounts")

# Remove constant rows
remove_constant_rows <- function(mat) { mat[rowSums(mat != 0) > 0, ]}

# logNorm_ref <- remove_constant_rows(assay(ref_sce, "logNorm"))
logNorm_sample <- remove_constant_rows(assay(sample_sce, "logNorm"))


# ## Dimension reduction
# pca_ref <- stats::prcomp(t(logNorm_ref), center = TRUE, scale. = TRUE) # reducedDim(ref_sce, "matPCs") <- matPCs(ref_sce, assay = "logNorm", 0.7)
pca_sample <- stats::prcomp(t(logNorm_sample), center = TRUE, scale. = TRUE) # reducedDim(sample_sce, "matPCs") <- matPCs(sample_sce, assay = "logNorm", 0.7)

# reducedDim(ref_sce, "matPCs") <- pca_ref$x
reducedDim(sample_sce, "matPCs") <- pca_sample$x

# Cell types
cellTypes <- sample_sce[["cellTypes"]]


# Check for NA values and remove them from ref_sce
valid_cells <- !is.na(cellTypes) & complete.cases(reducedDim(sample_sce, "matPCs"))
cellTypes <- cellTypes[valid_cells]
sample_sce <- sample_sce[, valid_cells]

# Check for balanced classes and adjust if necessary in sample_sce
cellType_counts <- table(cellTypes)
LimNum <- 5
if (any(cellType_counts < LimNum)) {
  warning(paste0("Some classes in reference have fewer than ",LimNum," samples. This may cause issues with the classifier."))
  min_class <- names(cellType_counts[cellType_counts < LimNum])
  for (cls in min_class) {
    sample_sce <- sample_sce[, cellTypes != cls]
    cellTypes <- cellTypes[cellTypes != cls]
  }
}


# Run scReClassify
set.seed(123)
cellTypes.reclassify <- multiAdaSampling(sample_sce, cellTypes, reducedDimName = "matPCs",
                                         classifier = Set_classifier, percent = Set_percent, L = Set_L)


# Save the new annotations to meta.data
Query_Seurat$label_scReClassify <- cellTypes.reclassify$final
