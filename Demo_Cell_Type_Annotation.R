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

rm(pbmc.rna, pbmc.atac)

# Set Actual_Cell_Type in seuratObject_Sample
seuratObject_Sample@meta.data[["Actual_Cell_Type"]] <- seuratObject_Sample@meta.data[["seurat_annotations"]]
seuratObject_Sample@meta.data$Actual_Cell_Type <- gsub("_", "  ", seuratObject_Sample@meta.data$Actual_Cell_Type)

# Set Actual_Cell_Type in seuratObject_Ref
seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]
seuratObject_Ref@meta.data$Actual_Cell_Type <- gsub("_", "  ", seuratObject_Ref@meta.data$Actual_Cell_Type)


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
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_singleR_NoReject", "label_singleR")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_singleR_NoReject')

## scmap
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_scmap_NoReject", "label_scmap")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scmap_NoReject')

## SCINA
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_SCINA_NoReject", "label_SCINA")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_SCINA_NoReject')

## scPred
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_scPred_NoReject", "label_scPred")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scPred_NoReject')

## CHETAH
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_CHETAH_NoReject", "label_CHETAH")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_CHETAH_NoReject')

## scClassify
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_scClassify_NoReject", "label_scClassify")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scClassify_NoReject')

## Seurat
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                             "label_Seurat_NoReject", "label_Seurat")
seuratObject_Sample <- FUN_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_Seurat_NoReject')

#### Visualization ####
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_singleR_NoReject")

DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_singleR_ConfStat")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_scmap_ConfStat")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_SCINA_ConfStat")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_scPred_ConfStat")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_CHETAH_ConfStat")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_scClassify_ConfStat")
DimPlot(seuratObject_Sample, reduction = "umap", group.by = "label_Seurat_ConfStat")

source("Set_plot_color.R")
source("PlotFun_Histogram.R")
metadata <- seuratObject_Sample@meta.data %>% as.data.frame()

## singleR
Plot_Hist_singleR_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_singleR_ConfStat', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)
Plot_Hist_singleR_Count

Plot_Hist_singleR_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_singleR_ConfStat', Note_Title = "",
                                   type = "proportion", color_vector = color_Class)
Plot_Hist_singleR_Prop

## scmap
Plot_Hist_scmap_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scmap_ConfStat', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)
Plot_Hist_scmap_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scmap_ConfStat', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## SCINA
Plot_Hist_SCINA_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_SCINA_ConfStat', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_SCINA_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_SCINA_ConfStat', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## scPred
Plot_Hist_scPred_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scPred_ConfStat', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_scPred_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scPred_ConfStat', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)


## CHETAH
Plot_Hist_CHETAH_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_CHETAH_ConfStat', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_CHETAH_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_CHETAH_ConfStat', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## scClassify
Plot_Hist_scClassify_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scClassify_ConfStat', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_scClassify_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_scClassify_ConfStat', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## Seurat
Plot_Hist_Seurat_Count <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_Seurat_ConfStat', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_Seurat_Prop <- plot_histogram(metadata, 'Actual_Cell_Type', 'label_Seurat_ConfStat', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## Sum
Plot_Hist_singleR_Count + Plot_Hist_scmap_Count + Plot_Hist_SCINA_Count +
Plot_Hist_scPred_Count + Plot_Hist_CHETAH_Count +Plot_Hist_scClassify_Count +
Plot_Hist_Seurat_Count

Plot_Hist_singleR_Prop + Plot_Hist_scmap_Prop + Plot_Hist_SCINA_Prop +
Plot_Hist_scPred_Prop + Plot_Hist_CHETAH_Prop +Plot_Hist_scClassify_Prop +
Plot_Hist_Seurat_Prop


################################################################################
#### VICTOR ####
# if(!require("devtools")) install.packages("devtools"); library(devtools)
# if(!require("VICTOR")) install_github("Charlene717/VICTOR"); library(VICTOR)

## singleR
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_singleR_NoReject")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_Diag(seuratObject_Sample, paste0("Diag_VICTOR_label_singleR_NoReject"),
                                                 paste0("DiagPara_VICTOR_label_singleR_NoReject"),
                                                 annotation_col = "label_singleR_NoReject")


## scmap
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_scmap_NoReject")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_Diag(seuratObject_Sample, paste0("Diag_VICTOR_label_scmap_NoReject"),
                                                 paste0("DiagPara_VICTOR_label_scmap_NoReject"),
                                                 annotation_col = "label_scmap_NoReject")

## SCINA
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_SCINA_NoReject")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_Diag(seuratObject_Sample, paste0("Diag_VICTOR_label_SCINA_NoReject"),
                                                 paste0("DiagPara_VICTOR_label_SCINA_NoReject"),
                                                 annotation_col = "label_SCINA_NoReject")

## scPred
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_scPred_NoReject")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_Diag(seuratObject_Sample, paste0("Diag_VICTOR_label_scPred_NoReject"),
                                                 paste0("DiagPara_VICTOR_label_scPred_NoReject"),
                                                 annotation_col = "label_scPred_NoReject")


## CHETAH
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_CHETAH_NoReject")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_Diag(seuratObject_Sample, paste0("Diag_VICTOR_label_CHETAH_NoReject"),
                                                 paste0("DiagPara_VICTOR_label_CHETAH_NoReject"),
                                                 annotation_col = "label_CHETAH_NoReject")

## scClassify
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_scClassify_NoReject")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_Diag(seuratObject_Sample, paste0("Diag_VICTOR_label_scClassify_NoReject"),
                                                 paste0("DiagPara_VICTOR_label_scClassify_NoReject"),
                                                 annotation_col = "label_scClassify_NoReject")


## Seurat
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_Seurat_NoReject")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_Diag(seuratObject_Sample, paste0("Diag_VICTOR_label_Seurat_NoReject"),
                                                 paste0("DiagPara_VICTOR_label_Seurat_NoReject"),
                                                 annotation_col = "label_Seurat_NoReject")


################################################################################
metadata <- seuratObject_Sample@meta.data %>% as.data.frame()

## singleR
Plot_Hist_singleR_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_singleR_NoReject', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)
Plot_Hist_singleR_Count_V

Plot_Hist_singleR_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_singleR_NoReject', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)
Plot_Hist_singleR_Prop_V

## scmap
Plot_Hist_scmap_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_scmap_NoReject', Note_Title = "",
                                        position_type = "stack", color_vector = color_Class)
Plot_Hist_scmap_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_scmap_NoReject', Note_Title = "",
                                       type = "proportion", color_vector = color_Class)

## SCINA
Plot_Hist_SCINA_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_SCINA_NoReject', Note_Title = "",
                                        position_type = "stack", color_vector = color_Class)

Plot_Hist_SCINA_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_SCINA_NoReject', Note_Title = "",
                                       type = "proportion", color_vector = color_Class)

## scPred
Plot_Hist_scPred_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_scPred_NoReject', Note_Title = "",
                                         position_type = "stack", color_vector = color_Class)

Plot_Hist_scPred_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_scPred_NoReject', Note_Title = "",
                                        type = "proportion", color_vector = color_Class)


## CHETAH
Plot_Hist_CHETAH_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_CHETAH_NoReject', Note_Title = "",
                                         position_type = "stack", color_vector = color_Class)

Plot_Hist_CHETAH_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_CHETAH_NoReject', Note_Title = "",
                                        type = "proportion", color_vector = color_Class)

## scClassify
Plot_Hist_scClassify_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_scClassify_NoReject', Note_Title = "",
                                             position_type = "stack", color_vector = color_Class)

Plot_Hist_scClassify_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_scClassify_NoReject', Note_Title = "",
                                            type = "proportion", color_vector = color_Class)

## Seurat
Plot_Hist_Seurat_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_Seurat_NoReject', Note_Title = "",
                                         position_type = "stack", color_vector = color_Class)

Plot_Hist_Seurat_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'DiagPara_VICTOR_label_Seurat_NoReject', Note_Title = "",
                                        type = "proportion", color_vector = color_Class)

## Sum
Plot_Hist_singleR_Count_V + Plot_Hist_scmap_Count_V + Plot_Hist_SCINA_Count_V +
  Plot_Hist_scPred_Count_V + Plot_Hist_CHETAH_Count_V +Plot_Hist_scClassify_Count_V +
  Plot_Hist_Seurat_Count_V

Plot_Hist_singleR_Prop_V + Plot_Hist_scmap_Prop_V + Plot_Hist_SCINA_Prop_V +
  Plot_Hist_scPred_Prop_V + Plot_Hist_CHETAH_Prop_V +Plot_Hist_scClassify_Prop_V +
  Plot_Hist_Seurat_Prop_V

################################################################################


Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)

# pdf(paste0(Name_ExportFolder,"/",Name_Export,"_",Set_AnnotM,"_",Set_ScoreM,"_AnnoDiagnosis_Hist.pdf"),
pdf(paste0(Name_time_wo_micro,"_AnnoDiagnosis_Hist.pdf"),
    width = 17, height = 17)

print(Plot_Hist_singleR_Count + Plot_Hist_scmap_Count + Plot_Hist_SCINA_Count +
        Plot_Hist_scPred_Count + Plot_Hist_CHETAH_Count +Plot_Hist_scClassify_Count +
        Plot_Hist_Seurat_Count)

print(Plot_Hist_singleR_Count_V + Plot_Hist_scmap_Count_V + Plot_Hist_SCINA_Count_V +
        Plot_Hist_scPred_Count_V + Plot_Hist_CHETAH_Count_V +Plot_Hist_scClassify_Count_V +
        Plot_Hist_Seurat_Count_V)

print(Plot_Hist_singleR_Prop + Plot_Hist_scmap_Prop + Plot_Hist_SCINA_Prop +
        Plot_Hist_scPred_Prop + Plot_Hist_CHETAH_Prop + Plot_Hist_scClassify_Prop +
        Plot_Hist_Seurat_Prop)

print(Plot_Hist_singleR_Prop_V + Plot_Hist_scmap_Prop_V + Plot_Hist_SCINA_Prop_V +
        Plot_Hist_scPred_Prop_V + Plot_Hist_CHETAH_Prop_V +Plot_Hist_scClassify_Prop_V +
        Plot_Hist_Seurat_Prop_V)

dev.off()



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

## -[] Function name
## -[] F1、Recall

## -[] VICTOR
## -[] VICTOR_glmnet
## -[V] VICTOR_StatROC
## -[] VICTOR_misc
## -[L] VICTOR_Prediction
## -[V] VICTOR_SeuratV5_SeuratV4

## -[] Fig Diag Accuracy

## -[] ATAC-seq
## -[] Sampling
## -[] Old sample
## -[] Big Data

## -[] UMAP

## -[] Reference Setting
## -[] Cell Type Name Function
## -[] Cell Type Name WO _
