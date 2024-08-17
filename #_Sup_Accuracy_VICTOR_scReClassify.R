##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("caret")) install.packages("caret"); library(caret)


if(!require("VICTOR")) devtools::install_github("Charlene717/VICTOR"); library(VICTOR)

if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("scReClassify")) BiocManager::install("scReClassify"); library(scReClassify)


#### Load Data ####
load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/#_Export_20240722/Export_GSE132044_20240720_scReClassify/20240712095055BUNQLI_MislabelB cell_Qry_10xV2_Ref_10xV2A/20240712095055BUNQLI.RData")







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


##### Visualization #####

#### Visualization CTAnnot label ####
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
  try({ DimPlot(seuratObject_Sample, reduction = "umap", group.by = label) })
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
plots <- lapply(labels_diag_para, function(label) {
  try({
    plot_histograms(metadata, 'Actual_Cell_Type', label, color_Class)
  })
})

# plots <- lapply(labels_diag_para, function(label) plot_histograms(metadata, 'Actual_Cell_Type', label, color_Class))

# 提取所有 Count 图
plots_count <- lapply(plots, `[[`, "Count")
gridExtra::grid.arrange(grobs = plots_count, ncol = 3)

# 提取所有 Prop 图
plots_prop <- lapply(plots, `[[`, "Prop")
gridExtra::grid.arrange(grobs = plots_prop, ncol = 3)


#### Visualization VICTOR ####
# 定義標籤
labels <- c(
  "label_singleR_NoReject",
  "label_scmap_NoReject",
  "label_SCINA_NoReject",
  "label_scPred_NoReject",
  "label_CHETAH_NoReject",
  "label_scClassify_NoReject",
  "label_Seurat_NoReject"
)

# 迭代绘图
plots_count_victor <- list()
plots_prop_victor <- list()
metadata <- seuratObject_Sample@meta.data %>% as.data.frame()

for (label in labels) {
  try({
    conf_stat_label <- paste0("ConfStat_VICTOR_", label)
    plots <- plot_histograms(metadata, "Actual_Cell_Type", conf_stat_label, color_Class)
    plots_count_victor <- c(plots_count_victor, list(plots[[1]]))
    plots_prop_victor <- c(plots_prop_victor, list(plots[[2]]))
  })
}

gridExtra::grid.arrange(grobs = plots_count_victor, ncol = 3)
gridExtra::grid.arrange(grobs = plots_prop_victor, ncol = 3)


#### Visualization scReClassify ####
# 迭代绘图
plots_count_scReClassify <- list()
plots_prop_scReClassify <- list()
metadata <- seuratObject_Sample@meta.data %>% as.data.frame()

for (label in labels) {
  try({
    conf_stat_label <- paste0("ReAnnot_scReClassify_", label,"_ConfStat")
    plots <- plot_histograms(metadata, "Actual_Cell_Type", conf_stat_label, color_Class)
    plots_count_scReClassify <- c(plots_count_scReClassify, list(plots[[1]]))
    plots_prop_scReClassify <- c(plots_prop_scReClassify, list(plots[[2]]))
  })
}

gridExtra::grid.arrange(grobs = plots_count_scReClassify, ncol = 3)
gridExtra::grid.arrange(grobs = plots_prop_scReClassify, ncol = 3)



# #### Export ####
# # Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
#
# ## Export PDF
# try({
#   pdf(paste0(export_directory,"/",Name_ExportFolder,"/",Name_Export,"_AnnoDiagnosis_Hist.pdf"),
#       # pdf(paste0(Name_time_wo_micro,"_AnnoDiagnosis_Hist.pdf"),
#       width = 17, height = 17)
#
#   # 绘制并输出图像
#   gridExtra::grid.arrange(grobs = plots_count, ncol = 3) %>% print()
#   gridExtra::grid.arrange(grobs = plots_count_victor, ncol = 3) %>% print()
#   gridExtra::grid.arrange(grobs = plots_count_scReClassify, ncol = 3) %>% print()
#
#   gridExtra::grid.arrange(grobs = plots_prop, ncol = 3) %>% print()
#   gridExtra::grid.arrange(grobs = plots_prop_victor, ncol = 3) %>% print()
#   gridExtra::grid.arrange(grobs = plots_prop_scReClassify, ncol = 3) %>% print()
#
#   dev.off()
#
# })
#
# ## Export MetaData
# write.table(data.frame(ID=rownames(seuratObject_Sample@meta.data), seuratObject_Sample@meta.data),
#             file=paste0(export_directory,"/",Name_ExportFolder,"/",Name_Export,"_metadataSamp.tsv"),
#             # file=paste0(Name_time_wo_micro,"_metadataSamp.tsv"),
#             quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')
# write.table(data.frame(ID=rownames(seuratObject_Ref@meta.data), seuratObject_Ref@meta.data),
#             file=paste0(export_directory,"/",Name_ExportFolder,"/",Name_Export,"_metadataRef.tsv"),
#             # file=paste0(Name_time_wo_micro,"_metadataRef.tsv"),
#             quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')
#
# # Remove Plot Object
# plot_objs <- grep("^[Pp]lot", ls(), value = TRUE)
# rm(list = plot_objs[sapply(plot_objs, function(obj) !is.function(get(obj)))])
#
# save.image(paste0(export_directory,"/",Name_ExportFolder,"/",Name_Export,".RData"))
#
