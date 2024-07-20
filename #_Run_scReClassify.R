# PMID31874628 scReClassify
# https://bioconductor.org/packages/release/bioc/vignettes/scReClassify/inst/doc/scReClassify.html


# ##### Presetting ######
# rm(list = ls()) # Clean variable ##* Comment out if Run All
# memory.limit(150000)


# load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_20240712/Export_GSE132044_MislabelB cell/20240712095055BUNQLI_MislabelB cell_Qry_10xV2_Ref_10xV2A/20240712095055BUNQLI.RData")

#### Load Packages ####
# Load necessary packages
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("caret")) install.packages("caret"); library(caret)

if (!require("scReClassify")) BiocManager::install("scReClassify"); library(scReClassify)
if (!require("DT")) install.packages("DT"); library(DT)
if (!require("mclust")) install.packages("mclust"); library(mclust)
if(!require("SingleCellExperiment"))  BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)


#### Run scReClassify ####
source("#_FUN_CellTypeAnnot.R")

# seuratObject_Sample <- Fun_scReClassify(seuratObject_Sample, Set_RefAnnoCol = "label_singleR", Set_classifier = "svm", Set_percent = 1, Set_L = 10)

# List of annotation columns to process
annoCols <- c("label_singleR_NoReject", "label_scmap_NoReject", "label_SCINA_NoReject",
              "label_scPred_NoReject", "label_CHETAH_NoReject", "label_scClassify_NoReject",
              "label_Seurat_NoReject")

# Apply run_scReClassify for each annotation column
for (annoCol in annoCols) {
  seuratObject_Sample <- Fun_scReClassify(seuratObject_Sample, Set_AnnoCol = annoCol, Set_classifier = "svm", Set_percent = 1, Set_L = 10)
}


#### DiagnosticMetrics ####
source("#_FUN_Metrics_CellTypeAnnot.R")

label_pairs <- list(
  c("label_singleR_NoReject", "ReAnnot_scReClassify_label_singleR_NoReject"),
  c("label_scmap_NoReject", "ReAnnot_scReClassify_label_scmap_NoReject"),
  c("label_SCINA_NoReject", "ReAnnot_scReClassify_label_SCINA_NoReject"),
  c("label_scPred_NoReject", "ReAnnot_scReClassify_label_scPred_NoReject"),
  c("label_CHETAH_NoReject", "ReAnnot_scReClassify_label_CHETAH_NoReject"),
  c("label_scClassify_NoReject", "ReAnnot_scReClassify_label_scClassify_NoReject"),
  c("label_Seurat_NoReject", "ReAnnot_scReClassify_label_Seurat_NoReject")
)


# 定義一個函數來處理每個標籤對
process_labels <- function(seuratObject, actual, no_reject, label) {
  seuratObject <- FUN_Confusion_Matrix(seuratObject, actual, no_reject, label)
  seuratObject <- FUN_CTAnnot_Accuracy(seuratObject, actual, no_reject)
  return(seuratObject)
}

# 依次處理每個標籤對
for (labels in label_pairs) {
  try({ seuratObject_Sample <- process_labels(seuratObject_Sample, "Actual_Cell_Type", labels[1], labels[2]) })
}
#
#
# #### Visualization ####
# # 定義標籤名稱列表
# labels_diag_para <- c(
#   "label_singleR_ConfStat",
#   "label_scmap_ConfStat",
#   "label_SCINA_ConfStat",
#   "label_scPred_ConfStat",
#   "label_CHETAH_ConfStat",
#   "label_scClassify_ConfStat",
#   "label_Seurat_ConfStat"
# )
#
# # 繪製 UMAP 圖
# for (label in labels_diag_para) {
#   try({ DimPlot(seuratObject_Sample, reduction = "umap", group.by = label) })
# }
#
# source("Set_plot_color.R")
# source("PlotFun_Histogram.R")
# metadata <- seuratObject_Sample@meta.data %>% as.data.frame()
#
# # 定義一個函數來處理每個標籤對
# plot_histograms <- function(metadata, actual, label, color_vector) {
#   list(
#     Count = plot_histogram(metadata, actual, label, Note_Title = "", position_type = "stack", color_vector = color_vector),
#     Prop = plot_histogram(metadata, actual, label, Note_Title = "", type = "proportion", color_vector = color_vector)
#   )
# }
#
# # 使用 lapply 来依次处理每个标签对
# plots <- lapply(labels_diag_para, function(label) {
#   try({
#     plot_histograms(metadata, 'Actual_Cell_Type', label, color_Class)
#   })
# })
#
# # plots <- lapply(labels_diag_para, function(label) plot_histograms(metadata, 'Actual_Cell_Type', label, color_Class))
#
# # 提取所有 Count 图
# plots_count <- lapply(plots, `[[`, "Count")
# gridExtra::grid.arrange(grobs = plots_count, ncol = 3)
#
# # 提取所有 Prop 图
# plots_prop <- lapply(plots, `[[`, "Prop")
# gridExtra::grid.arrange(grobs = plots_prop, ncol = 3)


