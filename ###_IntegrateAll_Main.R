## To-do list
## -[T] 去除Unassign
## -[T] 檢查程式寫錯的地方
## -[T] 重整同平台與跨平台的設定
## -[T] 整體數據檢查表
## -[] 精簡程式碼

## -[] 其他呈現方法(特別針對細胞種類去看，目前是一個sample為單位去算)

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
Dataset <- "scRNAseqPanc" # "GSE132044"  # "scRNAseqPanc"  # "HLCA_core"
Figure_Note <- Dataset

## Set export
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))

Name_Export <- paste0("IntegrateAll_",Name_FileID,"_",Dataset)

Name_ExportFolder <- paste0("Export_",Name_Export)
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder


#### Set Loading ###
## 設定主目錄
# main_directory <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_20240710_S"

if(Dataset == "GSE132044"){
  main_directory <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_20240712"
}else if(Dataset == "scRNAseqPanc"){
  main_directory <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_scRNAseqPanc_20240711"
}else if(Dataset == "HLCA_core"){
  main_directory <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results_HLCA_core/Export_HLCA_core_20240808"
}



# 獲取所有子目錄
subdirectories <- list.dirs(main_directory, recursive = TRUE, full.names = TRUE)

# 初始化一個空的dataframe
combined_data <- data.frame()

# 遍歷所有子目錄
for (subdir in subdirectories) {
  try({
    # 獲取當前子目錄中所有_metadataSamp.tsv結尾的檔案
    tsv_files <- list.files(subdir, pattern = "_metadataSamp\\.tsv$", full.names = TRUE)

    # 讀取並合併這些檔案
    for (file in tsv_files) {
      file_data <- read_tsv(file)
      combined_data <- bind_rows(combined_data, file_data)
    }
  })

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

# if(Set_VICTORS_Score == "SVDLRglmnet_ROC"){
#   # Remove rows with 'Method' containing 'EnsembleE2' or 'scPred' (except 'scPred' and 'scPred_SVGLRglmnet')
#   long_data <- long_data %>%
#     filter(!(grepl("EnsembleE2", Method) | (grepl("scPred", Method) & !Method %in% c("scPred", "scPred_SVGLRglmnet_ROC"))))
#
#   long_data <- long_data %>%
#     filter(!grepl("^[^_]*_[^_]*$", Method), # Removes rows with methods containing exactly one underscore
#            Actual_Cell_Type != "Unassigned") # Removes rows where Actual_Cell_Type is "Unassigned"
#
#   # 修改 long_data$Method 列，将包含 SVGLRglmnet_ROC 的方法名改为 VICTORS
#   long_data <- long_data %>%
#     mutate(Method = str_replace(Method, "SVGLRglmnet_ROC", "VICTORS"))
#
# }else if(Set_VICTORS_Score == "scPred_ROC"){
#   long_data <- long_data %>%
#     filter(Method %in% c("singleR", "scmap", "SCINA", "scPred", "singleR_scPred_ROC", "scmap_scPred_ROC", "SCINA_scPred_ROC", "scPred_scPred_ROC"))
#   long_data <- long_data %>%
#     filter(Actual_Cell_Type != "Unassigned") # Removes rows where Actual_Cell_Type is "Unassigned"
#   long_data <- long_data %>%
#     mutate(Method = str_replace(Method, "scPred_ROC", "VICTORS"))
#
# }else if(Set_VICTORS_Score =="EnsembleE2_ROC"){
#   long_data <- long_data %>%
#     filter(Method %in% c("singleR", "scmap", "SCINA", "scPred", "singleR_EnsembleE2_ROC", "scmap_EnsembleE2_ROC", "SCINA_EnsembleE2_ROC", "scPred_EnsembleE2_ROC"))
#   long_data <- long_data %>%
#     filter(Actual_Cell_Type != "Unassigned") # Removes rows where Actual_Cell_Type is "Unassigned"
#   long_data <- long_data %>%
#     mutate(Method = str_replace(Method, "EnsembleE2_ROC", "VICTORS"))
# }


## Remove duplicates
long_data <- long_data %>%
  dplyr::select(FileID, Sample_Platform, Ref_Platform, Mislabel_CellType,
                Actual_Cell_Type, Metric, Method,Value) %>%
  distinct()

long_data$Value <- as.numeric(long_data$Value)

# if (grepl("Baron", Set_load) && Set_Mislabel == "NoneMislabel") {
#   long_data <- long_data %>%
#     filter(!Actual_Cell_Type %in% c("Fibroblast cell","B cell","Endocrine cell"))
# }


# 绘制 boxplot
ggplot(long_data, aes(x = Mislabel_CellType, y = Value, fill = Method)) +
  geom_boxplot() +
  facet_wrap(~ Metric, scales = "free_y") + # 为不同的性能指标创建分面
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # 倾斜X轴标签以更好地显示
  labs(x = "Actual Cell Type", y = "Value", fill = "Annotation Method") +
  # scale_fill_brewer(palette = "Set1") # 使用色彩鲜明的调色板
  scale_fill_manual(values = color_Method) -> Plot_Box_All

Plot_Box_All

# plot_accuracy <- ggplot(long_data %>% filter(Metric == "Accuracy"),
#                         aes(x = Actual_Cell_Type, y = Value, fill = Method)) +
#   geom_boxplot() +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(x = "Actual Cell Type", y = "Accuracy", fill = "Annotation Method") +
#   scale_fill_brewer(palette = "Set1")
#
# plot_accuracy

custom_order <- c("singleR", "singleR_VICTOR", "scmap", "scmap_VICTOR",
                  "SCINA", "SCINA_VICTOR", "scPred","scPred_VICTOR",
                  "CHETAH", "CHETAH_VICTOR", "scClassify","scClassify_VICTOR",
                  "Seurat", "Seurat_VICTOR")
long_data$Method <- factor(long_data$Method, levels = custom_order)

custom_order2 <- long_data$Mislabel_CellType %>% unique() %>% sort()
custom_order2 <- custom_order2[custom_order2 !="None"]
custom_order2 <- c("None", custom_order2)
long_data$Mislabel_CellType <- factor(long_data$Mislabel_CellType, levels = custom_order2)



################################################################################
## Set Same_platform & Cross_platform
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)

# Create a function to group the same platforms
platform_group <- function(platform) {
  if (platform %in% c("10xV2", "10xV2A", "10xV2B")) {
    return("10xV2_Group")
  } else {
    return(platform)
  }
}

# Apply the function to long_data to add a column marking platform groups
long_data <- long_data %>%
  mutate(Sample_Platform_Group = sapply(Sample_Platform, platform_group),
         Ref_Platform_Group = sapply(Ref_Platform, platform_group))

# Same_platform
data_same_platform <- long_data %>%
  filter(Sample_Platform_Group == Ref_Platform_Group)

# Cross_platform
data_cross_platform <- long_data %>%
  filter(Sample_Platform_Group != Ref_Platform_Group)



## Boxplot
source("PlotFun_SetColor.R")
if (!exists("color_Method")) {color_Method <- setNames(character(0), character(0))}
color_Method <- update_color(c(data_same_platform$Method,"VICTOR"), color_Method)
if (!exists("color_legend")) {color_legend <- setNames(character(0), character(0))}
color_legend <- update_color(c(gsub(".*_VICTOR$", "VICTOR", as.character(data_same_platform$Method)),"VICTORS"), color_legend)

methods <- c("singleR", "scmap", "SCINA", "scPred", "CHETAH", "scClassify","Seurat" )
legend_set <- c(gsub(".*_VICTOR$", "VICTOR", as.character(long_data$Method))) %>% unique()




# Creating plots for different metrics and DataID
plot_accuracy_combined <- create_metric_plot(data_same_platform, "Accuracy", paste0(Figure_Note, " Accuracy by Missing Reference Cell Types - Same Platform"), color_Method , x_col = "Mislabel_CellType") /
  create_metric_plot(data_cross_platform, "Accuracy", paste0(Figure_Note, " Accuracy by Missing Reference Cell Types - Cross Platform"), color_Method, x_col = "Mislabel_CellType")
plot_accuracy_combined

plots_final_Accuracy_data_SamePlat <- create_and_combine_metric_plots(data_same_platform, methods, Figure_Note, "Accuracy", "Same Platform", 2, legend_set, color_legend, set_x_col = "Mislabel_CellType", set_ylimits = c(0.13, 1))
plots_final_Accuracy_data_SamePlat
plots_final_Accuracy_data_CrossPlat <- create_and_combine_metric_plots(data_cross_platform, methods, Figure_Note, "Accuracy", "Cross Platform", 2, legend_set, color_legend, set_x_col = "Mislabel_CellType", set_ylimits = c(0.13, 1))
plots_final_Accuracy_data_CrossPlat


plot_Recall_combined <- create_metric_plot(data_same_platform, "Recall", paste0(Figure_Note, " Recall by Missing Reference Cell Types - Same Platform"), color_Method, x_col = "Mislabel_CellType") /
  create_metric_plot(data_cross_platform, "Recall",  paste0(Figure_Note, " Recall by Missing Reference Cell Types - Cross Platform"), color_Method, x_col = "Mislabel_CellType")

plots_final_Recall_data_SamePlat <- create_and_combine_metric_plots(data_same_platform, methods, Figure_Note, "Recall", "Same Platform", 2, legend_set, color_legend, set_x_col = "Mislabel_CellType", set_ylimits = c(0.13, 1))
plots_final_Recall_data_CrossPlat <- create_and_combine_metric_plots(data_cross_platform, methods, Figure_Note, "Recall", "Cross Platform", 2, legend_set, color_legend, set_x_col = "Mislabel_CellType", set_ylimits = c(0.13, 1))


plot_Specificity_combined <- create_metric_plot(data_same_platform, "Specificity", paste0(Figure_Note, " Specificity by Missing Reference Cell Types - Same Platform"), color_Method, x_col = "Mislabel_CellType") /
  create_metric_plot(data_cross_platform, "Specificity",  paste0(Figure_Note, " Specificity by Missing Reference Cell Types - Cross Platform"), color_Method, x_col = "Mislabel_CellType")

plots_final_Specificity_data_SamePlat <- create_and_combine_metric_plots(data_same_platform, methods, Figure_Note, "Specificity", "Same Platform", 2, legend_set, color_legend, set_x_col = "Mislabel_CellType", set_ylimits = c(0.13, 1))
plots_final_Specificity_data_CrossPlat <- create_and_combine_metric_plots(data_cross_platform, methods, Figure_Note, "Specificity", "Cross Platform", 2, legend_set, color_legend, set_x_col = "Mislabel_CellType", set_ylimits = c(0.13, 1))



pdf(paste0(Name_ExportFolder, "/", Name_Export,"_MainResult.pdf"),
    width = 15, height = 15)  #GSE132044# width = 22, height = 14) #HLCA_core#  width = 30, height = 17) #scRNAseqPanc# width = 22, height = 13)

print(plots_final_Accuracy_data_SamePlat)
print(plots_final_Accuracy_data_CrossPlat)

print(plots_final_Recall_data_SamePlat)
print(plots_final_Recall_data_CrossPlat)


print(plots_final_Specificity_data_SamePlat)
print(plots_final_Specificity_data_CrossPlat)

dev.off()



# Set_Mislabel <- ""

pdf(paste0(Name_ExportFolder, "/", Name_Export,"_MainResult_Com.pdf"),
    width = 20, height = 20)

print(plot_accuracy_combined)
print(plot_Recall_combined)
print(plot_Specificity_combined)

print(Plot_Box_All)
dev.off()


# Remove Plot Object
plot_objs <- grep("^[Pp]lot", ls(), value = TRUE)
rm(list = plot_objs[sapply(plot_objs, function(obj) !is.function(get(obj)))])

# Export RData
save.image(paste0(Name_ExportFolder,"/", Name_Export,".RData"))
# save.image(paste0(Name_time_wo_micro,"_IntegrateAll.RData"))


################################################################################
################################################################################
#### Summary table ####

if(!require("dplyr")) install.packages("dplyr"); library(dplyr)

create_boxplot_table <- function(data, metric, value_var) {
  data %>%
    filter(Metric == metric) %>%
    group_by(Mislabel_CellType, Method) %>%
    summarize(
      Min = min(!!sym(value_var), na.rm = TRUE),
      Q1 = quantile(!!sym(value_var), 0.25, na.rm = TRUE),
      Median = median(!!sym(value_var), na.rm = TRUE),
      Q3 = quantile(!!sym(value_var), 0.75, na.rm = TRUE),
      Max = max(!!sym(value_var), na.rm = TRUE),
      .groups = 'drop'
    )
}

# # 调用函数
# Same_platform_accuracy_Quartile <- create_boxplot_table(data_same_platform, "Accuracy", "Value")
#
# # 查看结果
# print(Same_platform_accuracy_Quartile)



# 生成同平台数据集的准确度表格
try(TenX_platform_accuracy_Quartile <- create_boxplot_table(data_10x, "Accuracy", "Value"))
try(Same_platform_accuracy_Quartile <- create_boxplot_table(data_same_platform, "Accuracy", "Value"))
try(Cross_platform_accuracy_Quartile <- create_boxplot_table(data_cross_platform, "Accuracy", "Value"))

# 保存 TenX_platform_accuracy_Quartile 为 TSV 文件
try(write_tsv(TenX_platform_accuracy_Quartile, paste0(Name_ExportFolder, "/", Name_Export, "_10xPlatform.tsv")))

# 保存 Same_platform_accuracy_Quartile 为 TSV 文件
try(write_tsv(Same_platform_accuracy_Quartile, paste0(Name_ExportFolder, "/", Name_Export, "_SamePlatform.tsv")))

# 保存 Cross_platform_accuracy_Quartile 为 TSV 文件
try(write_tsv(Cross_platform_accuracy_Quartile, paste0(Name_ExportFolder, "/", Name_Export, "_CrossPlatform.tsv")))


# 仅保留 Metric 为 "Accuracy" 的行
try({
  Same_platform_accuracy <- data_same_platform %>%
    filter(Metric == "Accuracy")
})

try({
  Cross_platform_accuracy <- data_cross_platform %>%
    filter(Metric == "Accuracy")
})



# 移除环境中的其他对象
rm(list=setdiff(ls(), c("Name_ExportFolder", "Name_Export",
                        "Same_platform_accuracy", "Cross_platform_accuracy",
                        "Same_platform_accuracy_Quartile","Cross_platform_accuracy_Quartile")))

# Export RData
save.image(paste0(Name_ExportFolder,"/", Name_Export,"_S.RData"))


################################################################################
################################################################################
#### Data frame for check ####
selected_data_S <- selected_data %>%
  dplyr::select(FileID, Sample_Platform, Ref_Platform, Mislabel_CellType,
                contains("_Accuracy")) %>%
  distinct()

summary(is.na(selected_data_S))

