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

  # 调试信息：查看 metrics 数据框
  print(metrics)

  # 将计算后的指标合并回原始数据框
  data_with_metrics <- left_join(data, metrics, by = group_vars)

  return(data_with_metrics)
}

# ## Test function
# combined_dataTTT <- add_performance_metrics(combined_data,
#                                             c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
#                                             "label_singleR_ConfStat")


## singleR
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "label_singleR_ConfStat")

## scmap
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "label_scmap_ConfStat")

## SCINA
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "label_SCINA_ConfStat")
## scPred
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "label_scPred_ConfStat")
## CHETAH
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "label_CHETAH_ConfStat")
## scClassify
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "label_scClassify_ConfStat")

## Seurat
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "label_Seurat_ConfStat")

## VICTOR
## singleR
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "Diag_VICTOR_label_singleR_NoReject")

## scmap
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "Diag_VICTOR_label_scmap_NoReject")

## SCINA
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "Diag_VICTOR_label_SCINA_NoReject")
## scPred
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "Diag_VICTOR_label_scPred_NoReject")
## CHETAH
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "Diag_VICTOR_label_CHETAH_NoReject")
## scClassify
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "Diag_VICTOR_label_scClassify_NoReject")

## Seurat
combined_data <- add_performance_metrics(combined_data,
                                         c("FileID"), #c("FileID", "Mislabel_CellType", "Actual_Cell_Type"),
                                         "Diag_VICTOR_label_Seurat_NoReject")



## 同平台

## 不同平台


#### Visualization ####



