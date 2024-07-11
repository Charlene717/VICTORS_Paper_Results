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

if(!require("Rsamtools")) BiocManager::install("Rsamtools"); library(Rsamtools)
if(!require("Signac")) install.packages("Signac"); library(Signac)


#### Load Data ####
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")
#
# # seuratObject_Sample <- pbmc.rna
# # seuratObject_Ref <- pbmc.rna
# # seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]
#
#
# # Randomly split cells into two halves
# set.seed(123) # for reproducibility
# cells <- colnames(pbmc.rna)
# # sample_cells <- sample(cells, length(cells) / 2)
# # ref_cells <- setdiff(cells, sample_cells)
#
# sample_cells <- sample(cells, length(cells) / 5)
# ref_cells <- sample(cells, length(cells) / 5)
#
# # Create Seurat objects
# ## Qry
# if (!require("Signac")) install.packages("Signac")
# if (!require("EnsDb.Hsapiens.v86")) BiocManager::install("EnsDb.Hsapiens.v86")
# library(Signac)
# library(EnsDb.Hsapiens.v86)
#
# # 假设 pbmc.atac 是您的 ATAC Seurat 对象
# # 提取 ATAC 计数数据并获取基因组注释
# atac_counts <- GetAssayData(pbmc.atac, slot = "counts", assay = "ATAC")
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotations) <- "UCSC"
# genome(annotations) <- "hg38"
#
# # 找到 ATAC 峰值的基因组位置对应的基因符号
# peaks <- StringToGRanges(rownames(atac_counts), sep = c("-", "-"))
# annotation_overlap <- findOverlaps(peaks, annotations)
#
# # 创建基因符号列表并移除 NA 值的行
# gene_names <- rep(NA, length(peaks))
# gene_names[queryHits(annotation_overlap)] <- annotations$gene_name[subjectHits(annotation_overlap)]
# valid_genes <- !is.na(gene_names)
# atac_counts <- atac_counts[valid_genes, ]
# rownames(atac_counts) <- make.unique(gene_names[valid_genes])
#
# # 创建一个新的 Seurat 对象，指定 assay 为 "RNA"，并保留原有的 Metadata
# seuratObject_Sample <- CreateSeuratObject(counts = atac_counts, assay = "RNA")
# seuratObject_Sample@meta.data <- pbmc.atac@meta.data
#
#
# ## Ref
# seuratObject_Ref <- subset(pbmc.rna, cells = ref_cells)

##
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing_Integ.RData")

## Qry
seuratObject_Sample <- pbmc.integ.atac

## Ref
seuratObject_Ref <- pbmc.integ.rna

rm(pbmc.rna, pbmc.atac, pbmc.integ.rna, pbmc.integ.atac)

# Set Actual_Cell_Type in seuratObject_Sample
seuratObject_Sample@meta.data[["Actual_Cell_Type"]] <- seuratObject_Sample@meta.data[["seurat_annotations"]]
seuratObject_Sample@meta.data$Actual_Cell_Type <- gsub("_", "  ", seuratObject_Sample@meta.data$Actual_Cell_Type)

names(seuratObject_Sample@assays[["RNA"]]@layers)[names(seuratObject_Sample@assays[["RNA"]]@layers) == "data.ATAC"] <- "data"
# seuratObject_Sample@assays[["RNA"]]@layers[["counts"]] <- seuratObject_Sample@assays[["RNA"]]@layers[["data"]]
seuratObject_Sample@active.ident <- factor(seuratObject_Sample@active.ident, levels = "ATAC")
levels(seuratObject_Sample@active.ident) <- "RNA"


seuratObject_Sample@assays[["ATAC"]] <- NULL
seuratObject_Sample@assays[["ACTIVITY"]]  <- NULL

# 修改 @cells 中的 dimnames
dimnames(seuratObject_Sample@assays[["RNA"]]@cells) <- list(dimnames(seuratObject_Sample@assays[["RNA"]]@cells)[[1]], c("data.RNA", "scale.data.RNA"))
# 修改 @features 中的 dimnames
dimnames(seuratObject_Sample@assays[["RNA"]]@features) <- list(dimnames(seuratObject_Sample@assays[["RNA"]]@features)[[1]], c("data.RNA", "scale.data.RNA"))






# # 将 layers 中的每一项移动到 assays[["RNA"]] 下
# seuratObject_Sample@assays[["RNA"]]@counts <- seuratObject_Sample@assays[["RNA"]]@layers[["counts"]]
# seuratObject_Sample@assays[["RNA"]]@data <- seuratObject_Sample@assays[["RNA"]]@layers[["data"]]
# seuratObject_Sample@assays[["RNA"]]@scale.data <- seuratObject_Sample@assays[["RNA"]]@layers[["scale.data"]]
# seuratObject_Sample@assays[["RNA"]]@layers <- NULL

DimPlot(seuratObject_Sample, group.by = c("orig.ident", "seurat_annotations"))

# Set Actual_Cell_Type in seuratObject_Ref
seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]
seuratObject_Ref@meta.data$Actual_Cell_Type <- gsub("_", "  ", seuratObject_Ref@meta.data$Actual_Cell_Type)

names(seuratObject_Ref@assays[["RNA"]]@layers)[names(seuratObject_Ref@assays[["RNA"]]@layers) == "counts.RNA"] <- "counts"
names(seuratObject_Ref@assays[["RNA"]]@layers)[names(seuratObject_Ref@assays[["RNA"]]@layers) == "data.RNA"] <- "data"

DimPlot(seuratObject_Ref, group.by = c("orig.ident", "seurat_annotations"))

# ##### Data Preprocessing #####
# Seurat_Prepocessing <- function(seurat_obj, Num_PCA = 50 ,Set_nfeatures = 2000 ) {
#
#   seurat_obj <- seurat_obj  %>%
#     NormalizeData() %>%
#     FindVariableFeatures(nfeatures = Set_nfeatures) %>%
#     ScaleData() %>%
#     RunPCA(npcs = Num_PCA) %>%
#     FindNeighbors(dims = 1:Num_PCA) %>%
#     FindClusters() %>% # resolution = 0.5
#     # RunTSNE(dims = 1:Num_PCA) %>%
#     RunUMAP(dims = 1:Num_PCA)
# }
#
# seuratObject_Sample <- Seurat_Prepocessing(seuratObject_Sample)
#
# source("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Code_Rec/Run_harmony.R")

#### Run Cell Type Annotation ####
source("#_FUN_CellTypeAnnot.R")
DefaultAssay(seuratObject_Sample) <- "RNA"
DefaultAssay(seuratObject_Ref) <- "RNA"

## singleR
seuratObject_Sample <- Run_singleR(seuratObject_Sample, seuratObject_Ref,
                                   seurat_version = "V5M")

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
source("#_FUN_Metrics_CellTypeAnnot.R")

## singleR
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                            "label_singleR_NoReject", "label_singleR")
seuratObject_Sample <- FUN_CTAnnot_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_singleR_NoReject')

## scmap
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                            "label_scmap_NoReject", "label_scmap")
seuratObject_Sample <- FUN_CTAnnot_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scmap_NoReject')

## SCINA
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                            "label_SCINA_NoReject", "label_SCINA")
seuratObject_Sample <- FUN_CTAnnot_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_SCINA_NoReject')

## scPred
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                            "label_scPred_NoReject", "label_scPred")
seuratObject_Sample <- FUN_CTAnnot_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scPred_NoReject')

## CHETAH
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                            "label_CHETAH_NoReject", "label_CHETAH")
seuratObject_Sample <- FUN_CTAnnot_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_CHETAH_NoReject')

## scClassify
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                            "label_scClassify_NoReject", "label_scClassify")
seuratObject_Sample <- FUN_CTAnnot_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_scClassify_NoReject')

## Seurat
seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
                                            "label_Seurat_NoReject", "label_Seurat")
seuratObject_Sample <- FUN_CTAnnot_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_Seurat_NoReject')

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
                    AnnotCellTypeColumn = "label_singleR_NoReject",
                    seurat_version = "V5")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_DiagTools(seuratObject_Sample, paste0("Diag_VICTOR_label_singleR_NoReject"),
                                                      paste0("ConfStat_VICTOR_label_singleR_NoReject"),
                                                      annotation_col = "label_singleR_NoReject")


## scmap
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_scmap_NoReject",
                    seurat_version = "V5")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_DiagTools(seuratObject_Sample, paste0("Diag_VICTOR_label_scmap_NoReject"),
                                                      paste0("ConfStat_VICTOR_label_scmap_NoReject"),
                                                      annotation_col = "label_scmap_NoReject")

## SCINA
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_SCINA_NoReject",
                    seurat_version = "V5")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_DiagTools(seuratObject_Sample, paste0("Diag_VICTOR_label_SCINA_NoReject"),
                                                      paste0("ConfStat_VICTOR_label_SCINA_NoReject"),
                                                      annotation_col = "label_SCINA_NoReject")

## scPred
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_scPred_NoReject",
                    seurat_version = "V5")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_DiagTools(seuratObject_Sample, paste0("Diag_VICTOR_label_scPred_NoReject"),
                                                      paste0("ConfStat_VICTOR_label_scPred_NoReject"),
                                                      annotation_col = "label_scPred_NoReject")


## CHETAH
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_CHETAH_NoReject",
                    seurat_version = "V5")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_DiagTools(seuratObject_Sample, paste0("Diag_VICTOR_label_CHETAH_NoReject"),
                                                      paste0("ConfStat_VICTOR_label_CHETAH_NoReject"),
                                                      annotation_col = "label_CHETAH_NoReject")

## scClassify
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_scClassify_NoReject",
                    seurat_version = "V5")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_DiagTools(seuratObject_Sample, paste0("Diag_VICTOR_label_scClassify_NoReject"),
                                                      paste0("ConfStat_VICTOR_label_scClassify_NoReject"),
                                                      annotation_col = "label_scClassify_NoReject")


## Seurat
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_Seurat_NoReject",
                    seurat_version = "V5")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference

seuratObject_Sample <- FUN_Confusion_Matrix_DiagTools(seuratObject_Sample, paste0("Diag_VICTOR_label_Seurat_NoReject"),
                                                      paste0("ConfStat_VICTOR_label_Seurat_NoReject"),
                                                      annotation_col = "label_Seurat_NoReject")


################################################################################
metadata <- seuratObject_Sample@meta.data %>% as.data.frame()

## singleR
Plot_Hist_singleR_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_singleR_NoReject', Note_Title = "",
                                            position_type = "stack", color_vector = color_Class)
Plot_Hist_singleR_Count_V

Plot_Hist_singleR_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_singleR_NoReject', Note_Title = "",
                                           type = "proportion", color_vector = color_Class)
Plot_Hist_singleR_Prop_V

## scmap
Plot_Hist_scmap_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_scmap_NoReject', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)
Plot_Hist_scmap_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_scmap_NoReject', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## SCINA
Plot_Hist_SCINA_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_SCINA_NoReject', Note_Title = "",
                                          position_type = "stack", color_vector = color_Class)

Plot_Hist_SCINA_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_SCINA_NoReject', Note_Title = "",
                                         type = "proportion", color_vector = color_Class)

## scPred
Plot_Hist_scPred_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_scPred_NoReject', Note_Title = "",
                                           position_type = "stack", color_vector = color_Class)

Plot_Hist_scPred_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_scPred_NoReject', Note_Title = "",
                                          type = "proportion", color_vector = color_Class)


## CHETAH
Plot_Hist_CHETAH_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_CHETAH_NoReject', Note_Title = "",
                                           position_type = "stack", color_vector = color_Class)

Plot_Hist_CHETAH_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_CHETAH_NoReject', Note_Title = "",
                                          type = "proportion", color_vector = color_Class)

## scClassify
Plot_Hist_scClassify_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_scClassify_NoReject', Note_Title = "",
                                               position_type = "stack", color_vector = color_Class)

Plot_Hist_scClassify_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_scClassify_NoReject', Note_Title = "",
                                              type = "proportion", color_vector = color_Class)

## Seurat
Plot_Hist_Seurat_Count_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_Seurat_NoReject', Note_Title = "",
                                           position_type = "stack", color_vector = color_Class)

Plot_Hist_Seurat_Prop_V <- plot_histogram(metadata, 'Actual_Cell_Type', 'ConfStat_VICTOR_label_Seurat_NoReject', Note_Title = "",
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

