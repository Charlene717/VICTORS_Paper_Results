#### To-do list ####
# -[T] 修改命名
# -[T] 寫成函數
# -[T] 套用到ATAC
# -[ ] 整理程式碼

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require(cowplot)) install.packages("cowplot"); library(cowplot)

source("Plot_CellAnnot_UMAP_Box.R")


#### Load dataset ####
Dataset <- "GSE132044_MisLabelB"
# Dataset <- "SeuratATAC"

# load("D:/Dropbox/###_VUMC/##_Research/VICTORS/Figures/Figure1/Export_GSE132044_MislabelB cell/20231212125506KYHDNV_Multi/20231212125506KYHDNV.RData")
# export_folder <- "D:/Dropbox/###_VUMC/##_Research/VICTORS/Figures/Figure1/"
# export_name <- "VICTOR"

if(Dataset == "GSE132044_MisLabelB"){
  load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_20240712/Export_GSE132044_MislabelB cell/20240712095055BUNQLI_MislabelB cell_Qry_10xV2_Ref_10xV2A/20240712095055BUNQLI.RData")
  # load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_20240712/Export_GSE132044_MislabelB cell/20240712095816SIXEYV_MislabelB cell_Qry_10xV2_Ref_10xV3/20240712095816SIXEYV.RData")
  # load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_20240712/Export_GSE132044_MislabelB cell/20240712095420NCXIDQ_MislabelB cell_Qry_10xV2_Ref_10xV2B/20240712095420NCXIDQ.RData")

}else if(Dataset == "SeuratATAC"){
  load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_SeuratPBMC_ATAC/Export_SeuratPBMC_ATAC20240711192333.RData")

}


# source("FUN_Plot_Beautify_UMAP_Box.R")
# source("Set_plot_color.R")

DimPlot(seuratObject_Sample, reduction = "umap", group.by = "Actual_Cell_Type")
DimPlot(seuratObject_Ref, reduction = "umap", group.by = "Actual_Cell_Type")


#### Set Parameter ####

## Set export
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 10)
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))

Name_Export <- paste0("Fig1_",Name_FileID,"_",Dataset)

Name_ExportFolder <- paste0("Export_",Name_Export)
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder



#### Data pre-processing ####
Metadata <- seuratObject_Sample@meta.data

if ("FileID" %in% colnames(Metadata)) {
  Metadata <- Metadata %>%
    dplyr::select("FileID","Actual_Cell_Type","seurat_clusters",
                  "Sample_Platform","Ref_Platform", contains("ConfStat"))
} else {
  Metadata <- Metadata %>%
    dplyr::select("Actual_Cell_Type", #"Sample_Platform","Ref_Platform", # "seurat_clusters",
                  contains("ConfStat"))
}



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

new_colnames <- sapply(colnames(Metadata), convert_method_names)
colnames(Metadata) <- new_colnames

unique(long_metadata$Method) # 查看修改后的Method列




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
plot_accuracy_single <- function(accuracy_data, method, dataset, color_CellType = NULL) {
  plot <- ggplot(accuracy_data[accuracy_data$Method == method, ], aes(x = Actual_Cell_Type, y = Accuracy, fill = Actual_Cell_Type)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +  # 设置透明度
    geom_text(aes(label = sprintf("%.2f", Accuracy)), vjust = 0.5, hjust = 1, size = 6, angle = 90) +
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

  # 根据 color_CellType_Ref 是否存在来决定填充颜色
  if (!is.null(color_CellType)) {
    plot <- plot + scale_fill_manual(values = color_CellType)
  } else {
    # 使用灰色填充，确保为每个类都提供灰色
    unique_types <- unique(accuracy_data$Actual_Cell_Type)
    gray_palette <- rep(alpha("#9e9b9b", 0.8), length(unique_types))
    names(gray_palette) <- unique_types
    plot <- plot + scale_fill_manual(values = gray_palette)  # 使用灰色透明
  }

  return(plot)
}


# 計算所有方法的準確率並繪製圖表
# Method1 <- c("singleR", "scmap", "SCINA", "scPred")
Method1 <- methods[1:(length(methods)/2)]

All_accuracy_data1 <- do.call(rbind, lapply(Method1, function(alg) {
  calculate_accuracy(Metadata, alg, alg)
}))

# 为每个算法生成图表对象
plots1 <- lapply(Method1, function(alg) plot_accuracy_single(All_accuracy_data1, alg, Dataset))
# plots1 <- lapply(Method1, function(alg) plot_accuracy_single(All_accuracy_data1, alg, Dataset, color_CellType = color_CellType_Ref))

# 使用 grid.arrange 并排排列所有图表，设置 nrow = 1 以排成一行
grid.arrange(grobs = plots1, nrow = 1)
grid.arrange(grobs = c(plots[1:(length(plots)/2)],plots1), nrow = 2)

##
# Method2_VICTOR <- c("singleR_VICTOR", "scmap_VICTOR", "SCINA_VICTOR", "scPred_VICTOR")
Method2_VICTOR <- methods[(length(methods)/2 +1 ):length(methods)]

All_accuracy_data2_VICTOR <- do.call(rbind, lapply(Method2_VICTOR, function(alg) {
  calculate_accuracy(Metadata, alg, alg)
}))

# 为每个算法生成图表对象
plots2 <- lapply(Method2_VICTOR, function(alg) plot_accuracy_single(All_accuracy_data2_VICTOR, alg, Dataset))

# 使用 grid.arrange 并排排列所有图表，设置 nrow = 1 以排成一行
grid.arrange(grobs = plots2, nrow = 1)
grid.arrange(grobs = c(plots[(length(plots)/2+1):length(plots)],plots2), nrow = 2)





#### Export ####
## Export PDF
pdf_width <- if (grepl("^GSE132044", Dataset)) 28 else 40
pdf_height <- if (grepl("^GSE132044", Dataset)) 12 else 14

pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MainResult.pdf"),
    width = pdf_width, height = pdf_height)
grid.arrange(grobs = c(plots[1:(length(plots)/2)],plots1), nrow = 2)
grid.arrange(grobs = c(plots[(length(plots)/2+1):length(plots)],plots2), nrow = 2)
dev.off()


# pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MainResult_Fig1_Accuracy.pdf"),
#     width = 16, height = 8)
# print(grid.arrange(grobs = plots1, nrow = 1))
# print(grid.arrange(grobs = plots2, nrow = 1))
# dev.off()

## Export tsv
All_accuracy_data <- rbind(All_accuracy_data1, All_accuracy_data2_VICTOR)
try(write_tsv(All_accuracy_data, paste0(Name_ExportFolder, "/", Name_Export, "_Accuracy_Table.tsv")))


# # 移除环境中的其他对象
# rm(list=setdiff(ls(), c("All_accuracy_data","Metadata","Name_ExportFolder","Name_Export",
#                         "long_metadata", "seuratObject_Sample", "seuratObject_Ref")))
#
# # Export RData
# save.image(paste0(Name_ExportFolder,"/", Name_Export,".RData"))
#
#

################################################################################
library(ggplot2)
library(cowplot)
library(dplyr)
library(rlang)

# 提取所有子图的 x 轴标签（细胞类型）
cell_types <- unique(unlist(lapply(plots1, function(plot) {
  x_var <- as_label(plot$mapping$x)  # 使用 as_label 提取变量名
  plot$data[[x_var]]
})))

# 创建编号与细胞类型的对照表
cell_type_mapping <- data.frame(
  Number = seq_along(cell_types),
  Cell_Type = cell_types
)

# 定義一個函數來處理每個 plot 的更新
update_plot <- function(plot, index) {
  plot <- plot + scale_x_discrete(labels = setNames(cell_type_mapping$Number, cell_types)) +
    theme(axis.text.x = element_text(angle = 0, size = 24, face = "plain"),
          axis.text.y = element_text(size = 24, face = "plain"),
          axis.title = element_text(size = 24),
          text = element_text(size = 18, face = "bold"),
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +  # 完全移除边距，减少空隙
    labs(title = gsub(" on Actual_Cell_Type", "", plot$labels$title)) +
    theme(plot.title = element_text(size = 24, face = "bold"))  # 放大标题字体

  if (index != 1) {
    plot <- plot + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }

  plot <- plot + theme(axis.title.x = element_blank())

  return(plot)
}

# 更新 plots 的函数
update_plots <- function(plots) {
  lapply(seq_along(plots), function(i) update_plot(plots[[i]], i))
}

# 将 plots 分为两组并更新
half_length <- length(plots) / 2
plots1_1 <- update_plots(plots[1:half_length])
plots1 <- update_plots(plots1)

plots2_1 <- update_plots(plots[1:half_length])
plots2 <- update_plots(plots2)

# 合并并打印 combined plots
combine_and_print <- function(plots, ncol = 7) {
  combined_plot <- plot_grid(
    plotlist = plots,
    ncol = ncol,
    align = 'hv',
    axis = "tb",
    rel_widths = rep(1, length(plots)),
    rel_heights = rep(1, length(plots))
  )
  print(combined_plot)
  return(combined_plot)
}

combined_plots1_1 <- combine_and_print(plots1_1)
combined_plots2_1 <- combine_and_print(plots2_1)
combined_Percent_plot <- combine_and_print(plots1)
combined_Percent_plot2 <- combine_and_print(plots2)

# 打印组合的 plots
print(combined_plots1_1 / combined_Percent_plot)
print(combined_plots2_1 / combined_Percent_plot2)


# 输出编号与细胞类型的对照表
print(cell_type_mapping)



#### Export ####
## Export PDF
pdf_width <- if (grepl("^GSE132044", Dataset)) 28 else 40
pdf_height <- if (grepl("^GSE132044", Dataset)) 12 else 14

pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MainResult_2.pdf"),
    width = pdf_width, height = pdf_height)
print(combined_plots1_1 / combined_Percent_plot)
print(combined_plots2_1 / combined_Percent_plot2)
dev.off()



# # 移除环境中的其他对象
# rm(list=setdiff(ls(), c("All_accuracy_data","Metadata","Name_ExportFolder","Name_Export",
#                         "long_metadata", "seuratObject_Sample", "seuratObject_Ref")))
#
# # Export RData
# save.image(paste0(Name_ExportFolder,"/", Name_Export,".RData"))
#

