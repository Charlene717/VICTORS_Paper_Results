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


########################################

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
P1 <- do.call(gridExtra::grid.arrange, count_plots)
P1

# 提取並合併所有 Prop 圖
prop_plots <- lapply(plots, `[[`, "Prop")
P2 <- do.call(gridExtra::grid.arrange, prop_plots)
P2


