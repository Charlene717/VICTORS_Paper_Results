## To-do list
## -[ ] 折線圖
## -[ ] 精簡程式碼

## -[] 其他呈現方法

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

#### Set parameter ####
Dataset <- "GSE132044_B"
Figure_Note <- Dataset

## Set export
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))

Name_Export <- paste0("Sampling_",Name_FileID,"_",Dataset)

Name_ExportFolder <- paste0("Export_IntegrateAll_",Name_Export)
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder


#### Set Loading ###
## 設定主目錄

if(Dataset == "GSE132044_B"){
  main_directory <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_SamplingB_Seed123_20240711"
}else if(Dataset == ""){
  main_directory <- ""
}



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
# combined_data <- process_performance_metrics(combined_data, c("FileID", "Mislabel_CellType", "Actual_Cell_Type"), target_vars)

# head(combined_data)

## Remove duplicates
selected_data <- combined_data %>%
  dplyr::select(FileID, Sample_Platform, Ref_Platform, Mislabel_CellType,
                Actual_Cell_Type, contains("ConfStat")) %>%
  distinct()


################################################################################
#### Visualization ####
source("###_IntegAll_Visualization.R")

## long_data for all metrics # 将数据重塑为长格式，包括所有的指标
library(dplyr)
library(tidyr)
library(stringr)


long_data <- selected_data %>%
  pivot_longer(
    cols = ends_with(c("Accuracy", "Precision", "Recall", "F1", "Specificity")),
    names_to = "Metric_Method",
    values_to = "Value"
  ) %>%
  # Extract 'Metric' using regex
  mutate(Metric = str_extract(Metric_Method, "Accuracy|Precision|Recall|F1|Specificity")) %>%
  # Extract 'Method' using regex and remove unnecessary parts
  mutate(Method = str_remove(Metric_Method, "_?(Accuracy|Precision|Recall|F1|Specificity)$"),
         Method = str_remove(Method, "^ConfStat_VICTOR_"),
         Method = str_remove(Method, "^label_"),
         Method = str_replace(Method, "_NoReject", ""),
         Method = str_replace(Method, "_ConfStat", ""),
         Method = str_replace_all(Method, "_+", "_"),
         Method = str_remove(Method, "_$"),
         Method = if_else(str_detect(Metric_Method, "^ConfStat_VICTOR_"), paste0(Method, "_VICTOR"), Method)) %>%
  # Reformat 'Metric_Method' to 'Metric_Method'
  mutate(Metric_Method = paste0(Metric, "_", Method)) %>%
  # Convert 'Value' to numeric
  mutate(Value = as.numeric(Value)) %>%
  # Filter out rows where Metric is NA
  filter(!is.na(Metric)) %>%
  # Drop the original 'Metric_Method' column
  dplyr::select(-Metric_Method)

head(long_data) # Display the transformed data


## Remove duplicates
long_data <- long_data %>%
  dplyr::select(FileID, Sample_Platform, Ref_Platform, Mislabel_CellType,
                Actual_Cell_Type, Metric, Method,Value) %>%
  distinct()

long_data$Value <- as.numeric(long_data$Value)



# 绘制 折線圖
library(ggplot2)
library(dplyr)

# 过滤出Accuracy的数据
accuracy_data <- long_data %>%
  filter(Metric == "Accuracy")


# 定义颜色和线型
color_Method <- c(
  "singleR" = "#2862bf",
  "singleR_VICTOR" = "#2862bf",

  "scPred" = "#368a5b",
  "scPred_VICTOR" = "#368a5b",

  "SCINA" = "#e0385d",
  "SCINA_VICTOR" = "#e0385d",

  "scmap" = "#ad772f",
  "scmap_VICTOR" = "#ad772f",

  "CHETAH" = "#7ccf7f",
  "CHETAH_VICTOR" = "#7ccf7f",

  "scClassify" = "#8d389c",
  "scClassify_VICTOR" = "#8d389c",

  "Seurat" = "#42dbd9",
  "Seurat_VICTOR" = "#42dbd9"
)

linetype_Method <- c(
  "singleR" = "dashed",
  "singleR_VICTOR" = "solid",

  "scPred" = "dashed",
  "scPred_VICTOR" = "solid",

  "SCINA" = "dashed",
  "SCINA_VICTOR" = "solid",

  "scmap" = "dashed",
  "scmap_VICTOR" = "solid",

  "CHETAH" = "dashed",
  "CHETAH_VICTOR" = "solid",

  "scClassify" = "dashed",
  "scClassify_VICTOR" = "solid",

  "Seurat" = "dashed",
  "Seurat_VICTOR" = "solid"
)

# 按顺序排列Ref_Platform
accuracy_data <- long_data %>%
  filter(Metric == "Accuracy") %>%
  mutate(Ref_Platform = factor(Ref_Platform, levels = c("10xV2A_B0", "10xV2A_B10", "10xV2A_B20", "10xV2A_B30",
                                                        "10xV2A_B40", "10xV2A_B50", "10xV2A_B60", "10xV2A_B70",
                                                        "10xV2A_B80", "10xV2A_B90", "10xV2A_B100", "10xV2A_B110",
                                                        "10xV2A_B120", "10xV2A_B130", "10xV2A_B140")))

# 创建折线图
Plot_line <- ggplot(accuracy_data, aes(x = Ref_Platform, y = Value, color = Method, linetype = Method, group = Method)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = color_Method) +
  scale_linetype_manual(values = linetype_Method) +
  labs(title = "Accuracy across different Ref_Platforms by Method",
       x = "Ref Platform",
       y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Plot_line

#### Export ####
pdf(paste0(Name_ExportFolder, "/", Name_Export,"_MainResult.pdf"),
    width = 20, height = 20)

print(Plot_line)

dev.off()



# Export RData
save.image(paste0(Name_ExportFolder,"/", Name_Export,".RData"))
# save.image(paste0(Name_time_wo_micro,"_IntegrateAll.RData"))



