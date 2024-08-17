## Ref: https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette

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

#### Set Parameter ####

Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)

Name_ExportFolder <- "Export_SeuratPBMC_ATAC"
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder

Name_Export <- paste0(Name_ExportFolder,Name_time_wo_micro)



#### Load Data ####
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

seuratObject_Sample@meta.data[["nCount_RNA"]] <- NULL
seuratObject_Sample@meta.data[["nFeature_RNA"]] <- NULL
names(seuratObject_Sample@meta.data) <- gsub("_ATAC", "_RNA", names(seuratObject_Sample@meta.data))

# 修改 @cells 中的 dimnames
dimnames(seuratObject_Sample@assays[["RNA"]]@cells) <- list(dimnames(seuratObject_Sample@assays[["RNA"]]@cells)[[1]], c("data.RNA", "scale.data.RNA"))
# 修改 @features 中的 dimnames
dimnames(seuratObject_Sample@assays[["RNA"]]@features) <- list(dimnames(seuratObject_Sample@assays[["RNA"]]@features)[[1]], c("data.RNA", "scale.data.RNA"))


# 確保 data 層有正確的維度名稱
dimnames(seuratObject_Sample@assays[["RNA"]]@layers[["data"]]) <- list(
  rownames(seuratObject_Sample@assays[["RNA"]]),
  colnames(seuratObject_Sample@assays[["RNA"]])
)

# 確保 scale.data 層有正確的維度名稱
dimnames(seuratObject_Sample@assays[["RNA"]]@layers[["scale.data"]]) <- list(
  rownames(seuratObject_Sample@assays[["RNA"]]),
  colnames(seuratObject_Sample@assays[["RNA"]])
)

# Error # 将 layers 中的每一项移动到 assays[["RNA"]] 下
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
seuratObject_Sample <- Run_SCINA(seuratObject_Sample, seuratObject_Ref,
                                 seurat_version = "V5M")

## scPred
seuratObject_Sample <- Run_scPred(seuratObject_Sample, seuratObject_Ref)

## CHETAH
seuratObject_Sample <- Run_CHETAH(seuratObject_Sample, seuratObject_Ref)

## scClassify
seuratObject_Sample <- Run_scClassify(seuratObject_Sample, seuratObject_Ref,
                                      seurat_version = "V5M")

## Seurat
# seuratObject_Sample <- Run_Seurat_Annot(seuratObject_Sample, seuratObject_Ref)

# Rename the predicted.id column to label_Seurat_NoReject
colnames(seuratObject_Sample@meta.data)[which(colnames(seuratObject_Sample@meta.data) == "predicted.id")] <- "label_Seurat_NoReject"
seuratObject_Sample@meta.data$label_Seurat_NoReject <- gsub("_", "  ", seuratObject_Sample@meta.data$label_Seurat_NoReject)

# # Calculate the maximum prediction score for each cell
# seuratObject_Sample$mapping.score <- apply(seuratObject_Sample@meta.data[, grep("prediction.score", colnames(seuratObject_Sample@meta.data))], 1, max)
seuratObject_Sample$mapping.score <- seuratObject_Sample$prediction.score.max

# Add label_Seurat column based on label_Seurat_NoReject but mark cells with mapping.score < 0.8 as Unassign
seuratObject_Sample$label_Seurat <- ifelse(seuratObject_Sample$mapping.score < 0.8, "Unassign", seuratObject_Sample$label_Seurat_NoReject)
colnames(seuratObject_Sample@meta.data)[which(colnames(seuratObject_Sample@meta.data) == "mapping.score")] <- "label_Seurat_Scores"



## Pending
# ## scReClassify
# cellTypes.reclassify <- Fun_scReClassify(seuratObject_Sample, seuratObject_Ref)


################################################################################
#### DiagnosticMetrics ####
source("#_FUN_Metrics_CellTypeAnnot.R")

# ## Example of singleR
# seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
#                                             "label_singleR_NoReject", "label_singleR")
# seuratObject_Sample <- FUN_CTAnnot_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_singleR_NoReject')


# 定義標籤名稱列表
label_pairs <- list(
  c("label_singleR_NoReject", "label_singleR"),
  c("label_scmap_NoReject", "label_scmap"),
  c("label_SCINA_NoReject", "label_SCINA"),
  c("label_scPred_NoReject", "label_scPred"),
  c("label_CHETAH_NoReject", "label_CHETAH"),
  c("label_scClassify_NoReject", "label_scClassify"),
  c("label_Seurat_NoReject", "label_Seurat")
)

# 定義一個函數來處理每個標籤對
process_labels <- function(seuratObject, actual, no_reject, label) {
  seuratObject <- FUN_Confusion_Matrix(seuratObject, actual, no_reject, label)
  seuratObject <- FUN_CTAnnot_Accuracy(seuratObject, actual, no_reject)
  return(seuratObject)
}

# 依次處理每個標籤對
for (labels in label_pairs) {
  seuratObject_Sample <- process_labels(seuratObject_Sample, "Actual_Cell_Type", labels[1], labels[2])
}



################################################################################
#### Visualization ####
# 定義標籤名稱列表
labels_diag_para <- c(
  "label_singleR_ConfStat",
  "label_scmap_ConfStat",
  "label_SCINA_ConfStat",
  "label_scPred_ConfStat",
  "label_CHETAH_ConfStat",
  "label_scClassify_ConfStat",
  "label_Seurat_ConfStat"
)

# 繪製 UMAP 圖
for (label in labels_diag_para) {
  DimPlot(seuratObject_Sample, reduction = "umap", group.by = label)
}

source("Set_plot_color.R")
source("PlotFun_Histogram.R")
metadata <- seuratObject_Sample@meta.data %>% as.data.frame()

# 定義一個函數來處理每個標籤對
plot_histograms <- function(metadata, actual, label, color_vector) {
  list(
    Count = plot_histogram(metadata, actual, label, Note_Title = "", position_type = "stack", color_vector = color_vector),
    Prop = plot_histogram(metadata, actual, label, Note_Title = "", type = "proportion", color_vector = color_vector)
  )
}

# 使用 lapply 来依次处理每个标签对
plots <- lapply(labels_diag_para, function(label) plot_histograms(metadata, 'Actual_Cell_Type', label, color_Class))

# 提取所有 Count 图
plots_count <- lapply(plots, `[[`, "Count")
gridExtra::grid.arrange(grobs = plots_count, ncol = 3)

# 提取所有 Prop 图
plots_prop <- lapply(plots, `[[`, "Prop")
gridExtra::grid.arrange(grobs = plots_prop, ncol = 3)


################################################################################
################################################################################

#### VICTOR ####
if(!require("devtools")) install.packages("devtools"); library(devtools)
if(!require("VICTOR")) install_github("Charlene717/VICTOR"); library(VICTOR)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)

# ## Example of singleR
# VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
#                     ActualCellTypeColumn = "Actual_Cell_Type",
#                     AnnotCellTypeColumn = "label_singleR_NoReject")
#
# seuratObject_Sample <- VICTOR.lt$Query
# seuratObject_Ref <- VICTOR.lt$Reference
#
# seuratObject_Sample <- FUN_Confusion_Matrix_DiagTools(seuratObject_Sample, paste0("Diag_VICTOR_label_singleR_NoReject"),
#                                                       paste0("ConfStat_VICTOR_label_singleR_NoReject"),
#                                                       annotation_col = "label_singleR_NoReject")


# 定義標籤和診斷參數
labels <- c(
  "label_singleR_NoReject",
  "label_scmap_NoReject",
  "label_SCINA_NoReject",
  "label_scPred_NoReject",
  "label_CHETAH_NoReject",
  "label_scClassify_NoReject",
  "label_Seurat_NoReject"
)

# VICTOR和診斷處理函數
run_victor_and_diagnose <- function(seurat_obj_sample, seurat_obj_ref, actual_col, annot_col) {
  VICTOR.lt <- VICTOR(seurat_obj_sample, seurat_obj_ref,
                      ActualCellTypeColumn = actual_col,
                      AnnotCellTypeColumn = annot_col,
                      seurat_version = "V5")

  seurat_obj_sample <- VICTOR.lt$Query
  seurat_obj_ref <- VICTOR.lt$Reference

  diag_label <- paste0("Diag_VICTOR_", annot_col)
  conf_stat_label <- paste0("ConfStat_VICTOR_", annot_col)

  seurat_obj_sample <- FUN_Confusion_Matrix_DiagTools(seurat_obj_sample, diag_label, conf_stat_label, annotation_col = annot_col)

  return(list(seurat_obj_sample, seurat_obj_ref))
}

# 迭代執行VICTOR和診斷處理
for (label in labels) {
  result <- run_victor_and_diagnose(seuratObject_Sample, seuratObject_Ref, "Actual_Cell_Type", label)
  seuratObject_Sample <- result[[1]]
  seuratObject_Ref <- result[[2]]
}

# 迭代绘图
plots_count_victor <- list()
plots_prop_victor <- list()
metadata <- seuratObject_Sample@meta.data %>% as.data.frame()

for (label in labels) {
  conf_stat_label <- paste0("ConfStat_VICTOR_", label)
  plots <- plot_histograms(metadata, "Actual_Cell_Type", conf_stat_label, color_Class)
  plots_count_victor <- c(plots_count_victor, list(plots[[1]]))
  plots_prop_victor <- c(plots_prop_victor, list(plots[[2]]))
}

gridExtra::grid.arrange(grobs = plots_count_victor, ncol = 3)
gridExtra::grid.arrange(grobs = plots_prop_victor, ncol = 3)


#### Export ####

## Export PDF
try({
  pdf(paste0(Name_ExportFolder,"/",Name_Export,"_AnnoDiagnosis_Hist.pdf"),
      # pdf(paste0(Name_time_wo_micro,"_AnnoDiagnosis_Hist.pdf"),
      width = 17, height = 17)

  # 绘制并输出图像
  gridExtra::grid.arrange(grobs = plots_count, ncol = 3)
  gridExtra::grid.arrange(grobs = plots_count_victor, ncol = 3)
  gridExtra::grid.arrange(grobs = plots_prop, ncol = 3)
  gridExtra::grid.arrange(grobs = plots_prop_victor, ncol = 3)

  dev.off()

})

## Export MetaData
write.table(data.frame(ID=rownames(seuratObject_Sample@meta.data), seuratObject_Sample@meta.data),
            file=paste0(Name_ExportFolder,"/",Name_Export,"_metadataSamp.tsv"),
            # file=paste0(Name_time_wo_micro,"_metadataSamp.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')
write.table(data.frame(ID=rownames(seuratObject_Ref@meta.data), seuratObject_Ref@meta.data),
            file=paste0(Name_ExportFolder,"/",Name_Export,"_metadataRef.tsv"),
            # file=paste0(Name_time_wo_micro,"_metadataRef.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')

# Remove Plot Object
plot_objs <- grep("^[Pp]lot", ls(), value = TRUE)
rm(list = plot_objs[sapply(plot_objs, function(obj) !is.function(get(obj)))])

save.image(paste0(Name_ExportFolder,"/",Name_Export,".RData"))
# save.image(paste0(Name_time_wo_micro,".RData"))

# # Save small RData
# save(seuratObject_Sample, seuratObject_Ref, ROC_Summarize.lt,
#      file = paste0(Name_ExportFolder, "/", Name_Export, "_S.RData"))


