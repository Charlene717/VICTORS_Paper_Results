#### To-do list ####
# -[ ] 整理程式碼
# -[ ] 寫成函數
# -[ ] 寫入主程式


################################################################################
# PMID31874628 scReClassify
# https://bioconductor.org/packages/release/bioc/vignettes/scReClassify/inst/doc/scReClassify.html

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("caret")) install.packages("caret"); library(caret)

# Load necessary packages
if (!require("scReClassify")) BiocManager::install("scReClassify"); library(scReClassify)
if (!require("DT")) install.packages("DT"); library(DT)
if (!require("mclust")) install.packages("mclust"); library(mclust)
if(!require("SingleCellExperiment"))  BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)

## Set para
Set_Sam_Delet_Unknown <- TRUE
Set_Ref_Delet_Unknown <- TRUE
seurat_version <- "V5"

#### Load data ####
## Load sample
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Processed_Sample_OriQC/GSE132044_10xV2_Sample_CellNum1681_Seed123_Processed.RData") #; rm(Name_Export_o,Name_ExportFolder_o)

if(is.null(seuratObject@meta.data$`Actual_Cell_Type`)){
  seuratObject@meta.data$`Actual_Cell_Type` <- seuratObject@meta.data$`Cell_Type`
}
seuratObject <- UpdateSeuratObject(seuratObject)
Query_Seurat <- seuratObject; rm(seuratObject)
if(Set_Sam_Delet_Unknown){
  Query_Seurat <- subset(Query_Seurat, subset = Actual_Cell_Type != "Unknown")
  # Query_Seurat <- subset(Query_Seurat, subset = Cell_Type != "Unknown")
}

## SetIdent for Query_Seurat
try({ Query_Seurat <- Query_Seurat %>% SetIdent(value = "Annotation") })
DimPlot(Query_Seurat, label = TRUE, repel = TRUE) + NoLegend()
# DimPlot(Query_Seurat, reduction = "umap", group.by = "Annotation")

## Load reference
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Processed_Ref_OriQC/GSE132044_10xV2A_Ref_CellNum1611_Seed123_Processed.RData") #; rm(Name_Export_o,Name_ExportFolder_o)
seuratObject <- UpdateSeuratObject(seuratObject)
Reference_Seurat <- seuratObject; rm(seuratObject)
if(Set_Ref_Delet_Unknown){ Reference_Seurat <- subset(Reference_Seurat, subset = Actual_Cell_Type != "Unknown") }


## SetIdent for Reference_Seurat
Reference_Seurat <- Reference_Seurat %>% SetIdent(value = "Actual_Cell_Type")
DimPlot(Reference_Seurat, label = TRUE, repel = TRUE) + NoLegend()
# DimPlot(Reference_Seurat, reduction = "umap", group.by = "Actual_Cell_Type")

#### scReClassify ####
# Convert Seurat objects to SingleCellExperiment
# ref_sce <- as.SingleCellExperiment(Reference_Seurat)
sample_sce <- as.SingleCellExperiment(Query_Seurat)

# Standardize the data and create 'logNorm' assay
# ref_sce <- scater::logNormCounts(ref_sce)
sample_sce <- scater::logNormCounts(sample_sce)

# # Create 'logNorm' layer in assays
# SummarizedExperiment::assay(ref_sce, "logNorm") <- SummarizedExperiment::assay(ref_sce, "logcounts")
SummarizedExperiment::assay(sample_sce, "logNorm") <- SummarizedExperiment::assay(sample_sce, "logcounts")

# Remove constant rows
remove_constant_rows <- function(mat) { mat[rowSums(mat != 0) > 0, ]}

# logNorm_ref <- remove_constant_rows(assay(ref_sce, "logNorm"))
logNorm_sample <- remove_constant_rows(SummarizedExperiment::assay(sample_sce, "logNorm"))


# ## Dimension reduction
# pca_ref <- stats::prcomp(t(logNorm_ref), center = TRUE, scale. = TRUE) # reducedDim(ref_sce, "matPCs") <- matPCs(ref_sce, assay = "logNorm", 0.7)
pca_sample <- stats::prcomp(t(logNorm_sample), center = TRUE, scale. = TRUE) # reducedDim(sample_sce, "matPCs") <- matPCs(sample_sce, assay = "logNorm", 0.7)

# reducedDim(ref_sce, "matPCs") <- pca_ref$x
reducedDim(sample_sce, "matPCs") <- pca_sample$x

## (!!) Set Cell types
sample_sce[["cellTypes"]] <- sample_sce[["Actual_Cell_Type"]]
names(sample_sce@colData@listData[["cellTypes"]]) <- row.names(Query_Seurat@meta.data)
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
Set_classifier = "svm"
Set_percent = 1
Set_L = 10
# trace("multiAdaSampling", edit=TRUE) # layer = "data"
cellTypes.reclassify <- multiAdaSampling(sample_sce, cellTypes, reducedDimName = "matPCs",
                                         classifier = Set_classifier, percent = Set_percent, L = Set_L)


# Save the new annotations to meta.data
Query_Seurat$ReAnnot_scReClassify <- cellTypes.reclassify$final
DimPlot(Query_Seurat, label = TRUE, repel = TRUE , group.by = "ReAnnot_scReClassify") + NoLegend()
DimPlot(Query_Seurat, label = TRUE, repel = TRUE , group.by = "Actual_Cell_Type") + NoLegend()

# #### ConfStat ####
# source("#_FUN_Metrics_CellTypeAnnot.R")
# FUN_Confusion_Matrix(Query_Seurat,"Actual_Cell_Type", "Annotation","ReAnnot_scReClassify")
#



