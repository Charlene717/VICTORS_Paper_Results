#### Run Cell Type Annotation ####
source("#_FUN_CellTypeAnnot.R")

# 定義標註函數列表
annotation_functions <- list(
  Run_singleR,
  Run_scmap,
  Run_SCINA,
  Run_scPred,
  Run_CHETAH,
  Run_scClassify,
  Run_Seurat_Annot
)

# 依次執行每個標註函數
for (func in annotation_functions) {
  seuratObject_Sample <- func(seuratObject_Sample, seuratObject_Ref)
}


# ## singleR
# seuratObject_Sample <- Run_singleR(seuratObject_Sample, seuratObject_Ref)

################################################################################
#### DiagnosticMetrics ####
source("FUN_Metrics_CellTypeAnnot.R")

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


# ## singleR
# seuratObject_Sample <- FUN_Confusion_Matrix(seuratObject_Sample, "Actual_Cell_Type",
#                                             "label_singleR_NoReject", "label_singleR")
# seuratObject_Sample <- FUN_CTAnnot_Accuracy(seuratObject_Sample, 'Actual_Cell_Type', 'label_singleR_NoReject')

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

# 使用 lapply 來依次處理每個標籤對
plots <- lapply(labels_diag_para, function(label) plot_histograms(metadata, 'Actual_Cell_Type', label, color_Class))

# 提取並合併所有 Count 圖
count_plots <- lapply(plots, `[[`, "Count")
Plots_His_count_CTAnnot <- do.call(gridExtra::grid.arrange, count_plots)
Plots_His_count_CTAnnot

# 提取並合併所有 Prop 圖
prop_plots <- lapply(plots, `[[`, "Prop")
Plots_His_prop_CTAnnot <- do.call(gridExtra::grid.arrange, prop_plots)
Plots_His_prop_CTAnnot


################################################################################
################################################################################

#### VICTOR ####
# if(!require("devtools")) install.packages("devtools"); library(devtools)
# if(!require("VICTOR")) install_github("Charlene717/VICTOR"); library(VICTOR)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)


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
                      AnnotCellTypeColumn = annot_col)

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

# 繪圖函數
plot_histograms <- function(metadata, actual_col, conf_stat_label, color_vector) {
  plot_hist_count <- plot_histogram(metadata, actual_col, conf_stat_label, Note_Title = "", position_type = "stack", color_vector = color_vector)
  plot_hist_prop <- plot_histogram(metadata, actual_col, conf_stat_label, Note_Title = "", type = "proportion", color_vector = color_vector)

  return(list(plot_hist_count, plot_hist_prop))
}

# 迭代繪圖
plots_count <- list()
plots_prop <- list()
metadata <- seuratObject_Sample@meta.data %>% as.data.frame()

for (label in labels) {
  conf_stat_label <- paste0("ConfStat_VICTOR_", label)
  plots <- plot_histograms(metadata, "Actual_Cell_Type", conf_stat_label, color_Class)
  plots_count <- c(plots_count, list(plots[[1]]))
  plots_prop <- c(plots_prop, list(plots[[2]]))
}

# 合併所有的計數和比例圖
do.call(gridExtra::grid.arrange, plots_count)
do.call(gridExtra::grid.arrange, plots_prop)





