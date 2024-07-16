#### To-do list ####
# -[ ] 修改命名
# -[ ]
# -[ ] 寫成函數
# -[ ] 套用到ATAC

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require(cowplot)) install.packages("cowplot"); library(cowplot)

source("Plot_CellAnnot_UMAP_Box.R")

#### Load dataset ####
# load("D:/Dropbox/###_VUMC/##_Research/VICTORS/Figures/Figure1/Export_GSE132044_MislabelB cell/20231212125506KYHDNV_Multi/20231212125506KYHDNV.RData")
# export_folder <- "D:/Dropbox/###_VUMC/##_Research/VICTORS/Figures/Figure1/"
# export_name <- "VICTOR"

load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_20240712/Export_GSE132044_MislabelB cell/20240712095055BUNQLI_MislabelB cell_Qry_10xV2_Ref_10xV2A/20240712095055BUNQLI.RData")
export_folder <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_Figure1/"
export_name <- "VICTOR"

# source("FUN_Plot_Beautify_UMAP_Box.R")
# source("Set_plot_color.R")

DimPlot(seuratObject_Sample, reduction = "umap", group.by = "Actual_Cell_Type")
DimPlot(seuratObject_Ref, reduction = "umap", group.by = "Actual_Cell_Type")


# seuratObject_Sample@meta.data$singleR_VICTOR <- seuratObject_Sample@meta.data$DiagPara_label_singleR_NoReject_SVGLRglmnet_ROC
# seuratObject_Sample@meta.data$scmap_VICTOR <- seuratObject_Sample@meta.data$DiagPara_label_scmap_NoReject_SVGLRglmnet_ROC
# seuratObject_Sample@meta.data$SCINA_VICTOR <- seuratObject_Sample@meta.data$DiagPara_label_SCINA_NoReject_SVGLRglmnet_ROC
# seuratObject_Sample@meta.data$scPred_VICTOR <- seuratObject_Sample@meta.data$DiagPara_label_scPred_NoReject_Annot_SVGLRglmnet_ROC
#
# seuratObject_Sample@meta.data$singleR <- seuratObject_Sample@meta.data$label_singleR_DiagPara
# seuratObject_Sample@meta.data$scmap <- seuratObject_Sample@meta.data$label_scmap_DiagPara
# seuratObject_Sample@meta.data$SCINA <- seuratObject_Sample@meta.data$label_SCINA_DiagPara
# seuratObject_Sample@meta.data$scPred <- seuratObject_Sample@meta.data$label_scPred_DiagPara_Annot
# seuratObject_Sample@meta.data$ID <- seuratObject_Sample@meta.data$NAME
#
# Metadata <- seuratObject_Sample@meta.data
# Metadata <- Metadata[,c("ID","Actual_Cell_Type","seurat_clusters",
#                         "singleR_VICTOR","scmap_VICTOR","SCINA_VICTOR","scPred_VICTOR",
#                         "singleR","scmap","SCINA","scPred")]

Metadata <- seuratObject_Sample@meta.data
Metadata <- Metadata %>%
  dplyr::select("FileID","Actual_Cell_Type","seurat_clusters",
                 "Sample_Platform","Ref_Platform", contains("ConfStat"))



#### Plot Confusion Matrix Proportions ####
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)
if(!require('tidyr')) install.packages('tidyr'); library(tidyr)

# 将Metadata数据转换为长格式，便于处理和绘图
long_metadata <- Metadata %>%
  pivot_longer(cols = contains("ConfStat"),
               names_to = "Method",
               values_to = "CM_Value") %>%
  mutate(CM_Value = factor(CM_Value, levels = c("TP", "TN", "FP", "FN", "Other")))


## 修改Method名稱
library(dplyr)

# 定义转换函数
convert_method_names <- function(method) {
  if (grepl("^label_", method)) {
    method <- gsub("^label_", "", method)
    method <- gsub("_ConfStat$", "", method)
  } else if (grepl("^ConfStat_VICTOR_label_", method)) {
    method <- gsub("^ConfStat_VICTOR_label_", "", method)
    method <- gsub("_NoReject$", "", method)
    method <- paste0(method, "_VICTOR")
  }
  return(method)
}

# 应用转换函数到long_metadata的Method列
long_metadata <- long_metadata %>%
  mutate(Method = sapply(Method, convert_method_names))

# 查看修改后的Method列
unique(long_metadata$Method)



color_Class <- c(
  "TP" = "#b58b2a",
  "TN" = "#e8bc56",
  "FN" = "#368a5b",
  "FP" = "#73bd94",
  "Other" = "#919191"
)



plot_cm_values <- function(data, method) {
  filtered_data <- data %>% filter(Method == method)

  ggplot(filtered_data, aes(x = Actual_Cell_Type, fill = CM_Value)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = color_Class) +
    theme_minimal() +
    theme(text = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 16),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1) +
    labs(x = "Cell Type", y = "Proportion", title = paste("", method))
}

# 方法列表
methods <- unique(long_metadata$Method)

# 为每个方法生成图表对象
plots <- lapply(methods, function(method) plot_cm_values(long_metadata, method))

# 使用grid.arrange并排排列所有图表
grid.arrange(grobs = plots, nrow = 2)


#### Plot Accuracy ####
# 載入需要的包
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)

# 修改計算準確率的函數以保留方法信息
calculate_accuracy <- function(data, prediction_column, method) {
  accuracy_data <- data %>%
    group_by(Actual_Cell_Type) %>%
    summarise(
      TP = sum(.data[[prediction_column]] == "TP", na.rm = TRUE),
      TN = sum(.data[[prediction_column]] == "TN", na.rm = TRUE),
      FP = sum(.data[[prediction_column]] == "FP", na.rm = TRUE),
      FN = sum(.data[[prediction_column]] == "FN", na.rm = TRUE),
      Accuracy = (TP + TN) / (TP + TN + FP + FN),
      .groups = 'drop'
    )
  accuracy_data$Method <- method
  return(accuracy_data)
}


color_CellType_Ref <- list(
  "CD4+ T cell" = "#1f77b4",  # blue
  "B cell" = "#ff7f0e",       # orange
  "CD14+ monocyte" = "#2ca02c",  # green
  "Natural killer cell" = "#d62728",  # red
  "Cytotoxic T cell" = "#9467bd",  # purple
  "Megakaryocyte" = "#8c564b",  # brown
  "CD16+ monocyte" = "#e377c2",  # pink
  "Unknown" = "#7f7f7f",  # grey
  "unknown" = "#7f7f7f",  # grey
  "None" ="#7f7f7f",
  "Unassigned" = "#84e8d4",
  "Unassign" = "#84e8d4",
  "pruned" = "#f0ad8d",
  "Dendritic cell" = "#bcbd22",  # olive
  "Plasmacytoid dendritic cell" = "#17becf"  # light blue
)

# 在绘图函数中设置 alpha 参数
plot_accuracy_single <- function(accuracy_data, method) {
  plot <- ggplot(accuracy_data[accuracy_data$Method == method, ], aes(x = Actual_Cell_Type, y = Accuracy, fill = Actual_Cell_Type)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +  # 设置透明度
    geom_text(aes(label = sprintf("%.3f", Accuracy)), vjust = 0.5, hjust = 1, size = 6, angle = 90) +
    scale_fill_manual(values = color_CellType_Ref) +
    theme_minimal() +
    theme(text = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 16),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1) +
    ylim(0, 1.05) +
    labs(x = "Cell Type", y = "Accuracy", title = method)
  return(plot)
}


# 計算所有方法的準確率並繪製圖表
# Method1 <- c("singleR", "scmap", "SCINA", "scPred")
Method1 <- methods[1:(length(methods)/2)]

All_accuracy_data1 <- do.call(rbind, lapply(Method1, function(alg) {
  calculate_accuracy(Metadata, alg, alg)
}))

# 为每个算法生成图表对象
plots1 <- lapply(Method1, function(alg) plot_accuracy_single(All_accuracy_data1, alg))

# 使用 grid.arrange 并排排列所有图表，设置 nrow = 1 以排成一行
grid.arrange(grobs = plots1, nrow = 1)
grid.arrange(grobs = c(plots[1:4],plots1), nrow = 2)

##
Method2_VICTOR <- c("singleR_VICTOR", "scmap_VICTOR", "SCINA_VICTOR", "scPred_VICTOR")
Method2_VICTOR <- methods[(length(methods)/2 +1 ):length(methods)]

All_accuracy_data2_VICTOR <- do.call(rbind, lapply(Method2_VICTOR, function(alg) {
  calculate_accuracy(Metadata, alg, alg)
}))

# 为每个算法生成图表对象
plots2 <- lapply(Method2_VICTOR, function(alg) plot_accuracy_single(All_accuracy_data2_VICTOR, alg))

# 使用 grid.arrange 并排排列所有图表，设置 nrow = 1 以排成一行
grid.arrange(grobs = plots2, nrow = 1)
grid.arrange(grobs = c(plots[5:8],plots2), nrow = 2)





#### Export ####
## Export PDF
pdf(paste0(export_folder, "/", export_name, "_MainResult_Fig1.pdf"),
    width = 16, height = 12)
grid.arrange(grobs = c(plots[1:4],plots1), nrow = 2)
grid.arrange(grobs = c(plots[5:8],plots2), nrow = 2)
dev.off()


# pdf(paste0(export_folder, "/", export_name, "_MainResult_Fig1_Accuracy.pdf"),
#     width = 16, height = 8)
# print(grid.arrange(grobs = plots1, nrow = 1))
# print(grid.arrange(grobs = plots2, nrow = 1))
# dev.off()

## Export tsv
All_accuracy_data <- rbind(All_accuracy_data1, All_accuracy_data2_VICTOR)
try(write_tsv(All_accuracy_data, paste0(export_folder, "/", export_name, "_Accuracy_Table.tsv")))


# # 移除环境中的其他对象
# rm(list=setdiff(ls(), c("All_accuracy_data","Metadata","export_folder","export_name")))
#
# # Export RData
# save.image(paste0(export_folder,"/", export_name,"_Fig1_Accuracy.RData"))
