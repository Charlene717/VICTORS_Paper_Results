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

## Remove duplicates
selected_data <- combined_data %>%
  dplyr::select(FileID, Sample_DataID, Ref_DataID, Sample_DataID, Ref_DataID,
                Actual_Cell_Type, contains("ConfStat")) %>%
  distinct()


################################################################################
#### Visualization ####
source("###_IntegAll_Visualization.R")

## long_data for all metrics # 将数据重塑为长格式，包括所有的指标
library(dplyr)
library(tidyr)
library(stringr)

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


long_data <- long_data %>%
  dplyr::select(FileID, Sample_DataID, Ref_DataID, Sample_DataID, Ref_DataID,
                Actual_Cell_Type, Metric, Method,Value) %>%
  distinct()

long_data$Value <- as.numeric(long_data$Value)

# if (grepl("Baron", Set_load) && Set_Mislabel == "NoneMislabel") {
#   long_data <- long_data %>%
#     filter(!Actual_Cell_Type %in% c("Fibroblast cell","B cell","Endocrine cell"))
# }


# 绘制 boxplot
ggplot(long_data, aes(x = Actual_Cell_Type, y = Value, fill = Method)) +
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

custom_order2 <- long_data$Actual_Cell_Type %>% unique() %>% sort()
custom_order2 <- custom_order2[custom_order2 !="None"]
custom_order2 <- c("None", custom_order2)
long_data$Actual_Cell_Type <- factor(long_data$Actual_Cell_Type, levels = custom_order2)



################################################################################
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
color_Method <- update_color(c(data_same_DataID$Method,"VICTOR"), color_Method)
if (!exists("color_legend")) {color_legend <- setNames(character(0), character(0))}
color_legend <- update_color(c(gsub(".*_VICTOR$", "VICTOR", as.character(data_same_DataID$Method)),"VICTORS"), color_legend)

methods <- c("singleR", "scmap", "SCINA", "scPred", "CHETAH", "scClassify","Seurat" )
legend_set <- c(gsub(".*_VICTOR$", "VICTOR", as.character(long_data$Method))) %>% unique()


Figure_Note <- ""

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



# pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MainResult.pdf"),
#     width = 10, height = 10)

Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)

pdf(paste0(Name_time_wo_micro,"_MainResult.pdf"),
    width = 10, height = 20)

print(plots_final_Accuracy_data_SamePlat)
print(plots_final_Accuracy_data_CrossPlat)

print(plots_final_Recall_data_SamePlat)
print(plots_final_Recall_data_CrossPlat)


print(plots_final_Specificity_data_SamePlat)
print(plots_final_Specificity_data_CrossPlat)

dev.off()



Set_Mislabel <- ""
# pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MainResult_Com.pdf"),
#     width = 17, height = 17) #  width = 17, height = 17)

pdf(paste0(Name_time_wo_micro,"_MainResult_Com.pdf"),
    width = 10, height = 20)

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
# save.image(paste0(Name_ExportFolder,"/", Name_Export,".RData"))
save.image(paste0(Name_time_wo_micro,"_IntegrateAll.RData"))
