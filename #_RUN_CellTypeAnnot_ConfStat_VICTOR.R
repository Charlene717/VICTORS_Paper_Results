### Parameter ###
source("##_RunAll_CTAEvaluator_Set_Parameter.R")

## Export parameter information to txt file
writeLines(readLines("##_RunAll_CTAEvaluator_Set_Parameter.R"),
           con = paste0(Name_ExportFolder,"/",Name_Export,"_Parameter_Settings_Record.txt"))

##### Load data* #####
## Load sample
load(Path_Sample) #; rm(Name_Export_o,Name_ExportFolder_o)

if(is.null(seuratObject@meta.data$`Actual_Cell_Type`)){
  seuratObject@meta.data$`Actual_Cell_Type` <- seuratObject@meta.data$`Cell_Type`
}

seuratObject_Sample <- seuratObject; rm(seuratObject)
if(Set_Sam_Delet_Unknown){
  seuratObject_Sample <- subset(seuratObject_Sample, subset = Actual_Cell_Type != "Unknown")
  # seuratObject_Sample <- subset(seuratObject_Sample, subset = Cell_Type != "Unknown")
}

## SetIdent for seuratObject_Sample
seuratObject_Sample <- seuratObject_Sample %>% SetIdent(value = "Annotation")
DimPlot(seuratObject_Sample, label = TRUE, repel = TRUE) + NoLegend()
# DimPlot(seuratObject_Sample, reduction = "umap", group.by = "Annotation")

## Load reference
load(Path_Ref) #; rm(Name_Export_o,Name_ExportFolder_o)
seuratObject_Ref <- seuratObject; rm(seuratObject)
if(Set_Ref_Delet_Unknown){ seuratObject_Ref <- subset(seuratObject_Ref, subset = Actual_Cell_Type != "Unknown") }


if(Set_Ref_Delet_CTMetric){
  seuratObject_Ref <- subset(seuratObject_Ref, subset = Actual_Cell_Type != Set_Ref_Delet)
}


## SetIdent for seuratObject_Ref
seuratObject_Ref <- seuratObject_Ref %>% SetIdent(value = "Actual_Cell_Type")
DimPlot(seuratObject_Ref, label = TRUE, repel = TRUE) + NoLegend()
# DimPlot(seuratObject_Ref, reduction = "umap", group.by = "Actual_Cell_Type")


##### Data Preprocessing #####
## Seurat object Prepocessing
source("FUN_Seurat_Prepocessing.R")
if(Set_Run_Process_Sample){ seuratObject_Sample <- Seurat_Prepocessing(seuratObject_Sample, Num_PCA = Num_PCA ,Set_nfeatures = Set_nfeatures) }
if(Set_Run_Process_Ref){ seuratObject_Ref <- Seurat_Prepocessing(seuratObject_Ref, Num_PCA = Num_PCA ,Set_nfeatures = Set_nfeatures) }



#### Run Cell Type Annotation ####
source("#_FUN_CellTypeAnnot.R")

# ## Example of singleR
# seuratObject_Sample <- Run_singleR(seuratObject_Sample, seuratObject_Ref)


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
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)

## Export PDF
# pdf(paste0(Name_ExportFolder,"/",Name_Export,"_",Set_AnnotM,"_",Set_ScoreM,"_AnnoDiagnosis_Hist.pdf"),
pdf(paste0(Name_time_wo_micro,"_AnnoDiagnosis_Hist.pdf"),
    width = 17, height = 17)

# 绘制并输出图像
gridExtra::grid.arrange(grobs = plots_count, ncol = 3)
gridExtra::grid.arrange(grobs = plots_prop, ncol = 3)
gridExtra::grid.arrange(grobs = plots_count_victor, ncol = 3)
gridExtra::grid.arrange(grobs = plots_prop_victor, ncol = 3)

dev.off()

## Export MetaData
write.table(data.frame(ID=rownames(seuratObject_Sample@meta.data), seuratObject_Sample@meta.data),
            # file=paste0(Name_ExportFolder,"/",Name_Export,"_metadataSamp.tsv"),
            file=paste0(Name_time_wo_micro,"_metadataSamp.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')
write.table(data.frame(ID=rownames(seuratObject_Ref@meta.data), seuratObject_Ref@meta.data),
            # file=paste0(Name_ExportFolder,"/",Name_Export,"_metadataRef.tsv"),
            file=paste0(Name_time_wo_micro,"_metadataRef.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')

# Remove Plot Object
plot_objs <- grep("^[Pp]lot", ls(), value = TRUE)
rm(list = plot_objs[sapply(plot_objs, function(obj) !is.function(get(obj)))])

# save.image(paste0(Name_ExportFolder,"/",Name_Export,".RData"))
save.image(paste0(Name_time_wo_micro,".RData"))

# # Save small RData
# save(seuratObject_Sample, seuratObject_Ref, ROC_Summarize.lt,
#      file = paste0(Name_ExportFolder, "/", Name_Export, "_S.RData"))


