##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

if(!require("caret")) install.packages("caret"); library(caret)
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("readr")) install.packages("readr"); library(readr)

#### Set Loading ###
# 設定主目錄
main_directory <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_20240710_S"

# 獲取所有子目錄
subdirectories <- list.dirs(main_directory, recursive = TRUE, full.names = TRUE)

# 初始化一個空的dataframe
combined_data <- data.frame()

# 遍歷所有子目錄
for (subdir in subdirectories) {
  # 獲取當前子目錄中所有_metadataSamp.tsv結尾的檔案
  tsv_files <- list.files(subdir, pattern = "_metadataSamp\\.tsv$", full.names = TRUE)

  # 讀取並合併這些檔案
  for (file in tsv_files) {
    file_data <- read_tsv(file)
    combined_data <- bind_rows(combined_data, file_data)
  }
}

# 顯示整合後的dataframe
print(combined_data)




#### Visualization ####

#### Add Metrics ####
add_performance_metrics <- function(data, group_vars, target_var) {
  # 定义性能指标的名称
  accuracy_name <- paste0(target_var, "_Accuracy")
  precision_name <- paste0(target_var, "_Precision")
  recall_name <- paste0(target_var, "_Recall")
  f1_name <- paste0(target_var, "_F1")
  specificity_name <- paste0(target_var, "_Specificity")

  # 计算性能指标
  metrics <- data %>%
    group_by(across(all_of(group_vars))) %>%
    summarize(
      TP = sum(.data[[target_var]] == "TP", na.rm = TRUE),
      TN = sum(.data[[target_var]] == "TN", na.rm = TRUE),
      FP = sum(.data[[target_var]] == "FP", na.rm = TRUE),
      FN = sum(.data[[target_var]] == "FN", na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      !!accuracy_name := (TP + TN) / (TP + TN + FP + FN),
      !!precision_name := if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
      !!recall_name := if_else(TP + FN > 0, TP / (TP + FN), NA_real_),
      !!f1_name := if_else(TP + FP > 0 & TP + FN > 0, 2 * (TP / (TP + FP)) * (TP / (TP + FN)) / ((TP / (TP + FP)) + (TP / (TP + FN))), NA_real_),
      !!specificity_name := if_else(TN + FP > 0, TN / (TN + FP), NA_real_)
    ) %>%
    dplyr::select(-TP, -TN, -FP, -FN)  # 移除临时列

  # print(metrics) # 调试信息：查看 metrics 数据框

  # 将计算后的指标合并回原始数据框
  data_with_metrics <- left_join(data, metrics, by = group_vars)

  return(data_with_metrics)
}

# ## Test function
# combined_dataTTT <- add_performance_metrics(combined_data,
#                                             c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
#                                             "label_singleR_ConfStat")


if(!require("dplyr")) install.packages("dplyr"); library(dplyr)

# 定义一个函数来批量处理所有目标变量
process_performance_metrics <- function(data, group_vars, target_vars) {
  for (target_var in target_vars) {
    data <- add_performance_metrics(data, group_vars, target_var)
  }
  return(data)
}

# 定义所有目标变量
target_vars <- c("label_singleR_ConfStat", "label_scmap_ConfStat", "label_SCINA_ConfStat",
                 "label_scPred_ConfStat", "label_CHETAH_ConfStat", "label_scClassify_ConfStat",
                 "label_Seurat_ConfStat",
                 "ConfStat_VICTOR_label_singleR_NoReject", "ConfStat_VICTOR_label_scmap_NoReject",
                 "ConfStat_VICTOR_label_SCINA_NoReject", "ConfStat_VICTOR_label_scPred_NoReject",
                 "ConfStat_VICTOR_label_CHETAH_NoReject", "ConfStat_VICTOR_label_scClassify_NoReject",
                 "ConfStat_VICTOR_label_Seurat_NoReject")

# 执行批量处理
combined_data <- process_performance_metrics(combined_data, c("FileID"), target_vars)

# head(combined_data)


## Within-platform
combined_data_SamePlatform <- combined_data[combined_data$Sample_Platform == combined_data$Ref_Platform, ]

## Cross-platform
combined_data_CrossPlatform <- combined_data[combined_data$Sample_Platform != combined_data$Ref_Platform, ]



################################################################################
#### Visualization ####
source("###_IntegAll_Visualization.R")

# 將long_data篩選出不同組合的dataframe，Same_DataID、 Cross_DataID
# Same_DataID
data_same_DataID <- long_data %>%
  filter((Sample_DataID == Ref_DataID))

# Cross_DataID
data_cross_platform <- long_data %>%
  filter(!(Sample_DataID == Ref_DataID))


## Boxplot
source("PlotFun_SetColor.R")
if (!exists("color_Method")) {color_Method <- setNames(character(0), character(0))}
color_Method <- update_color(c(data_same_DataID$Method,"VICTORS"), color_Method)
if (!exists("color_legend")) {color_legend <- setNames(character(0), character(0))}
color_legend <- update_color(c(gsub(".*_VICTORS$", "VICTORS", as.character(data_same_DataID$Method)),"VICTORS"), color_legend)

methods <- c("singleR", "scmap", "SCINA", "scPred")
legend_set <- c(gsub(".*_VICTORS$", "VICTORS", as.character(long_data$Method))) %>% unique()

# Creating plots for different metrics and DataID
plot_accuracy_combined <- create_metric_plot(data_same_DataID, "Accuracy", paste0(Figure_Note, " Accuracy Across Actual Cell Types - Same DataID"), color_Method) /
  create_metric_plot(data_cross_platform, "Accuracy", paste0(Figure_Note, " Accuracy Across Actual Cell Types - Cross DataID"), color_Method)

plots_final_Accuracy_data_SamePlat <- create_and_combine_metric_plots(data_same_DataID, methods, Figure_Note, "Accuracy", "Same DataID", 2, legend_set, color_legend)
plots_final_Accuracy_data_CrossPlat <- create_and_combine_metric_plots(data_cross_platform, methods, Figure_Note, "Accuracy", "Cross DataID", 2, legend_set, color_legend)


plot_Recall_combined <- create_metric_plot(data_same_DataID, "Recall", paste0(Figure_Note, " Recall Across Actual Cell Types - Same DataID"), color_Method) /
  create_metric_plot(data_cross_platform, "Recall",  paste0(Figure_Note, " Recall Across Actual Cell Types - Cross DataID"), color_Method)

plots_final_Recall_data_SamePlat <- create_and_combine_metric_plots(data_same_DataID, methods, Figure_Note, "Recall", "Same DataID", 2, legend_set, color_legend)
plots_final_Recall_data_CrossPlat <- create_and_combine_metric_plots(data_cross_platform, methods, Figure_Note, "Recall", "Cross DataID", 2, legend_set, color_legend)


plot_Specificity_combined <- create_metric_plot(data_same_DataID, "Specificity", paste0(Figure_Note, " Specificity Across Actual Cell Types - Same DataID"), color_Method) /
  create_metric_plot(data_cross_platform, "Specificity",  paste0(Figure_Note, " Specificity Across Actual Cell Types - Cross DataID"), color_Method)

plots_final_Specificity_data_SamePlat <- create_and_combine_metric_plots(data_same_DataID, methods, Figure_Note, "Specificity", "Same DataID", 2, legend_set, color_legend)
plots_final_Specificity_data_CrossPlat <- create_and_combine_metric_plots(data_cross_platform, methods, Figure_Note, "Specificity", "Cross DataID", 2, legend_set, color_legend)



pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MainResult.pdf"),
    width = 10, height = 10)

print(plots_final_Accuracy_data_SamePlat)
print(plots_final_Accuracy_data_CrossPlat)

print(plots_final_Recall_data_SamePlat)
print(plots_final_Recall_data_CrossPlat)


print(plots_final_Specificity_data_SamePlat)
print(plots_final_Specificity_data_CrossPlat)

dev.off()


pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MainResult_Com.pdf"),
    width = 17, height = 17) #  width = 17, height = 17)
print(plot_accuracy_combined)
if(Set_Mislabel == "NoneMislabel"){
  print(plot_Recall_combined)
}else{
  print(plot_Specificity_combined)
}
print(Plot_Box_All)
dev.off()


# Remove Plot Object
plot_objs <- grep("^[Pp]lot", ls(), value = TRUE)
rm(list = plot_objs[sapply(plot_objs, function(obj) !is.function(get(obj)))])

# Export RData
save.image(paste0(Name_ExportFolder,"/", Name_Export,".RData"))




