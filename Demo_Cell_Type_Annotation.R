##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("caret")) install.packages("caret"); library(caret)

#### Load Function ####
if(!require("scPred")) devtools::install_github("powellgenomicslab/scPred"); library(scPred)
trace("project_query", edit=TRUE) # new_data <- GetAssayData(new, "data")[shared_features, ] # new_data <- GetAssayData(new, layer = "data")[shared_features, ]

source("#_FUN_CellTypeAnnot.R")


#### Load Data ####
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")

# seuratObject_Sample <- pbmc.rna
# seuratObject_Ref <- pbmc.rna
# seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]


# Randomly split cells into two halves
set.seed(123) # for reproducibility
cells <- colnames(pbmc.rna)
# sample_cells <- sample(cells, length(cells) / 2)
# ref_cells <- setdiff(cells, sample_cells)

sample_cells <- sample(cells, length(cells) / 5)
ref_cells <- sample(cells, length(cells) / 5)

# Create Seurat objects
seuratObject_Sample <- subset(pbmc.rna, cells = sample_cells)
seuratObject_Ref <- subset(pbmc.rna, cells = ref_cells)

# Set Actual_Cell_Type in seuratObject_Sample
seuratObject_Sample@meta.data[["Actual_Cell_Type"]] <- seuratObject_Sample@meta.data[["seurat_annotations"]]


# Set Actual_Cell_Type in seuratObject_Ref
seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]


#### Run Cell Type Annotation ####
## singleR
seuratObject_Sample <- Run_singleR(seuratObject_Sample, seuratObject_Ref)

## scmap
seuratObject_Sample <- Run_scmap(seuratObject_Sample, seuratObject_Ref)

## SCINA
seuratObject_Sample <- Run_SCINA(seuratObject_Sample, seuratObject_Ref)

## scPred
seuratObject_Sample <- Run_scPred(seuratObject_Sample, seuratObject_Ref)

## CHETAH
seuratObject_Sample <- Run_CHETAH(seuratObject_Sample, seuratObject_Ref)

## scClassify
seuratObject_Sample <- Run_scClassify(seuratObject_Sample, seuratObject_Ref)

## Seurat
seuratObject_Sample <- Run_Seurat_Annot(seuratObject_Sample, seuratObject_Ref)


## Pending
# ## scReClassify
# cellTypes.reclassify <- Fun_scReClassify(seuratObject_Sample, seuratObject_Ref)


#### DiagnosticMetrics ####
source("FUN_Metrics_CellTypeAnnot.R")

## singleR
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_singleR_NoReject", "label_singleR")

## scmap
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_scmap_NoReject", "label_scmap")

## SCINA
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_SCINA_NoReject", "label_SCINA")
## scPred
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_scPred_NoReject", "label_scPred")

## CHETAH
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_CHETAH_NoReject", "label_CHETAH")

## scClassify
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_scClassify_NoReject", "label_scClassify")

## Seurat
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_Seurat_NoReject", "label_Seurat")


#### Visualization ####
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_singleR_NoReject")

DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_singleR_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_scmap_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_SCINA_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_scPred_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_CHETAH_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_scClassify_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_Seurat_DiagPara")
