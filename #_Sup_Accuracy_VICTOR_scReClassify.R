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



# 確保所需的包已安裝並加載
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if(!require("stringr")) install.packages("stringr"); library(stringr)

# 提取以 _ConfStat 結尾或 ConfStat_ 開頭的列名
confstat_cols <- colnames(seuratObject_Sample@meta.data) %>%
  str_subset("(_ConfStat$|^ConfStat_)")

# 計算每個方法的 Accuracy
accuracy_data <- seuratObject_Sample@meta.data %>%
  select(all_of(confstat_cols)) %>%
  summarise(across(everything(), ~{
    TP <- sum(. == "TP", na.rm = TRUE)
    TN <- sum(. == "TN", na.rm = TRUE)
    FP <- sum(. == "FP", na.rm = TRUE)
    FN <- sum(. == "FN", na.rm = TRUE)
    Accuracy <- (TP + TN) / (TP + TN + FP + FN)
    return(Accuracy)
  }, .names = "Accuracy_{col}")) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Accuracy") %>%
  drop_na()  # 忽略NA值

# 將列名轉化為更易讀的格式
accuracy_data$Method <- gsub("_ConfStat|ConfStat_", "", accuracy_data$Method)

# 繪製直方圖
ggplot(accuracy_data, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Accuracy across different methods",
       x = "Method", y = "Accuracy")


# 保留包含 VICTOR 或 scReClassify 的列
accuracy_data <- accuracy_data %>%
  filter(str_detect(Method, "VICTOR") | str_detect(Method, "scReClassify"))

# 修改 Method 名稱
accuracy_data <- accuracy_data %>%
  mutate(Method = case_when(
    str_detect(Method, "VICTOR") ~ str_replace(Method, "Accuracy_VICTOR_label_(.*?)_NoReject", "\\1_VICTOR"),
    str_detect(Method, "scReClassify") ~ str_replace(Method, "Accuracy_ReAnnot_scReClassify_label_(.*?)_NoReject", "\\1_scReClassify"),
    TRUE ~ Method
  ))

# 設定顏色
color_mapping <- c("VICTOR" = "#335769", "scReClassify" = "#55a6cf")
accuracy_data$Diagnostic_Tool <- ifelse(str_detect(accuracy_data$Method, "VICTOR"), "VICTOR", "scReClassify")

# 提取副標題信息
subtitle_text <- paste0(
  "Query: ", seuratObject_Sample@misc[["BasicInfo"]][["Platform"]],
  "; Reference: ", seuratObject_Sample@misc[["CTAnnot"]][["Ref_Platform"]]
)

# 确保数据框非空
if(nrow(accuracy_data) == 0) {
  stop("accuracy_data is empty. Please check the filtering step.")
}

# 绘制直方图
ggplot(accuracy_data, aes(x = Method, y = Accuracy, fill = Diagnostic_Tool)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.5) +  # 设置黑色边框和较小的粗细
  scale_fill_manual(values = color_mapping) +
  theme_minimal(base_size = 20) +  # 设置基础字体大小
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 22),  # 放大x轴标签字体
    axis.text.y = element_text(size = 22),  # 放大y轴标签字体
    axis.title.x = element_text(size = 24, face = "bold"),  # 放大x轴标题字体
    axis.title.y = element_text(size = 24, face = "bold"),  # 放大y轴标题字体
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),  # 放大标题字体并居中
    plot.subtitle = element_text(size = 22, face = "italic", hjust = 0.5),  # 放大副标题字体并居中
    legend.title = element_text(size = 22),  # 放大图例标题字体
    legend.text = element_text(size = 20),  # 放大图例内容字体
    panel.border = element_rect(color = "black", fill = NA, size = 1.5)  # 增加黑色粗框
  ) +
  labs(
    title = "Accuracy across different methods",
    subtitle = subtitle_text,  # 添加副标题
    x = "Method", y = "Accuracy"
  ) + scale_y_continuous(limits = c(0,1))


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
