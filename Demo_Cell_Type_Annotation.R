##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("caret")) install.packages("caret"); library(caret)

#### Load Function ####
if(!require("devtools")) install.packages("devtools"); library(devtools)
if(!require("scPred")) devtools::install_github("powellgenomicslab/scPred"); library(scPred)
trace("project_query", edit=TRUE) # layer = "data"
if(!require("VICTOR")) devtools::install_github("Charlene717/VICTOR"); library(VICTOR)

# new_data <- GetAssayData(new, "data")[shared_features, ] # new_data <- GetAssayData(new, layer = "data")[shared_features, ]

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
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_singleR_NoReject')

## scmap
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_scmap_NoReject", "label_scmap")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scmap_NoReject')

## SCINA
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_SCINA_NoReject", "label_SCINA")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_SCINA_NoReject')

## scPred
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_scPred_NoReject", "label_scPred")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scPred_NoReject')

## CHETAH
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_CHETAH_NoReject", "label_CHETAH")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_CHETAH_NoReject')

## scClassify
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_scClassify_NoReject", "label_scClassify")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scClassify_NoReject')

## Seurat
seuratObject_Sample <- FUN_DiagnosticMetrics(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_Seurat_NoReject", "label_Seurat")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_Seurat_NoReject')

#### Visualization ####
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_singleR_NoReject")

DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_singleR_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_scmap_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_SCINA_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_scPred_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_CHETAH_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_scClassify_DiagPara")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_Seurat_DiagPara")

source("Set_plot_color.R")
source("PlotFun_Histogram.R")
metadata <- seuratObject_Sample@meta.data %>% as.data.frame()

## singleR
Plot_Hist_singleR_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_singleR_DiagPara', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)
Plot_Hist_singleR_Count

Plot_Hist_singleR_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_singleR_DiagPara', Note_Title = "",
                                   type = "proportion", color_vector = color_Class)
Plot_Hist_singleR_Prop

## scmap
Plot_Hist_scmap_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scmap_DiagPara', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)
Plot_Hist_scmap_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scmap_DiagPara', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## SCINA
Plot_Hist_SCINA_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_SCINA_DiagPara', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_SCINA_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_SCINA_DiagPara', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## scPred
Plot_Hist_scPred_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scPred_DiagPara', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_scPred_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scPred_DiagPara', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)


## CHETAH
Plot_Hist_CHETAH_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_CHETAH_DiagPara', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_CHETAH_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_CHETAH_DiagPara', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## scClassify
Plot_Hist_scClassify_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scClassify_DiagPara', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_scClassify_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scClassify_DiagPara', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## Seurat
Plot_Hist_Seurat_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_Seurat_DiagPara', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_Seurat_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_Seurat_DiagPara', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## Sum
Plot_Hist_singleR_Count + Plot_Hist_scmap_Count + Plot_Hist_SCINA_Count +
Plot_Hist_scPred_Count + Plot_Hist_CHETAH_Count +Plot_Hist_scClassify_Count +
Plot_Hist_Seurat_Count

Plot_Hist_singleR_Prop + Plot_Hist_scmap_Prop + Plot_Hist_SCINA_Prop +
Plot_Hist_scPred_Prop + Plot_Hist_CHETAH_Prop +Plot_Hist_scClassify_Prop +
Plot_Hist_Seurat_Prop

Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)

# pdf(paste0(Name_ExportFolder,"/",Name_Export,"_",Set_AnnotM,"_",Set_ScoreM,"_AnnoDiagnosis_Hist.pdf"),
pdf(paste0(Name_time_wo_micro,"_AnnoDiagnosis_Hist.pdf"),
    width = 17, height = 17)

print(Plot_Hist_singleR_Count + Plot_Hist_scmap_Count + Plot_Hist_SCINA_Count +
        Plot_Hist_scPred_Count + Plot_Hist_CHETAH_Count +Plot_Hist_scClassify_Count +
        Plot_Hist_Seurat_Count)

print(Plot_Hist_singleR_Prop + Plot_Hist_scmap_Prop + Plot_Hist_SCINA_Prop +
        Plot_Hist_scPred_Prop + Plot_Hist_CHETAH_Prop + Plot_Hist_scClassify_Prop +
        Plot_Hist_Seurat_Prop)

dev.off()

################################################################################
#### VICTOR ####
# if(!require("devtools")) install.packages("devtools"); library(devtools)
# if(!require("VICTOR")) install_github("Charlene717/VICTOR"); library(VICTOR)

## singleR
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_singleR_NoReject")


## scmap

## SCINA

## scPred


## CHETAH

## scClassify

## Seurat


################################################################################

## singleR

## scmap

## SCINA

## scPred


## CHETAH

## scClassify

## Seurat

################################################################################

#### To-do list ####
## -[] Metric_Other

## -[] VICTOR
## -[] VICTOR_glmnet
## -[] VICTOR_StatROC
## -[] VICTOR_misc
## -[] VICTOR_Prediction
## -[V] VICTOR_SeuratV5_SeuratV4

## -[] Fig Diag Accuracy

## -[] ATAC-seq
## -[] Sampling
## -[] Old sample
## -[] Big Data

## -[] UMAP

## -[] Reference Setting
## -[] Cell Type Name Function
