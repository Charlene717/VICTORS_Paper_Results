##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


##### Load Package #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)


##### Load Data #####
load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/#_Export_20240717/Export_Fig1_2024071601VAI_GSE132044_MisLabelB_Ref10xV2A/Fig1_2024071601VAI_GSE132044_MisLabelB_Fig1_Accuracy.RData")

##### Data Processing #####

#### Cell type number ####
# 從 Seurat 物件提取資料、計算總和並添加總和列
Count_ActualCT_Query.df <- seuratObject_Sample@meta.data$Actual_Cell_Type %>%
  table() %>%
  as.data.frame() %>%
  setNames(c("Cell_Type", "Freq")) %>%
  rbind(data.frame(Cell_Type = "Total", Freq = sum(.$Freq)))

# Count_ActualCT_Ref.df <- seuratObject_Ref@meta.data$Actual_Cell_Type %>%
#   table() %>%
#   as.data.frame() %>%
#   setNames(c("Cell_Type", "Freq")) %>%
#   rbind(data.frame(Cell_Type = "Total", Freq = sum(.$Freq)))

#### ConfStat ####
# # 確保所需的包已安裝並載入
# if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}
#
# # 將感興趣的數據提取到新變數中
# metadata <- seuratObject_Sample@meta.data
#
# # 計算不同細胞類型被不同標註結果標記的數量
# result_table <- metadata %>%
#   group_by(Actual_Cell_Type, label_singleR_NoReject, label_singleR_ConfStat) %>%
#   summarise(Cell_Count = n()) %>%
#   arrange(Actual_Cell_Type, label_singleR_NoReject, label_singleR_ConfStat)
#
# # 查看結果表格
# print(result_table)


# 確保所需的包已安裝並載入
if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}

# 提取meta.data
metadata <- seuratObject_Sample@meta.data

# 創建一個函數來處理每個標註結果
summarize_labeling <- function(label_col, conf_col) {
  metadata %>%
    group_by(Actual_Cell_Type, !!sym(label_col), !!sym(conf_col)) %>%
    summarise(Cell_Count = n()) %>%
    arrange(Actual_Cell_Type, !!sym(label_col), !!sym(conf_col)) %>%
    rename(Label = !!sym(label_col), ConfStat = !!sym(conf_col))
}


# 將所有結果存儲在一個列表中
results_list <- list(
  singleR = summarize_labeling("label_singleR_NoReject", "label_singleR_ConfStat"),
  scmap = summarize_labeling("label_scmap_NoReject", "label_scmap_ConfStat"),
  SCINA = summarize_labeling("label_SCINA_NoReject", "label_SCINA_ConfStat"),
  scPred = summarize_labeling("label_scPred_NoReject", "label_scPred_ConfStat"),
  CHETAH = summarize_labeling("label_CHETAH_NoReject", "label_CHETAH_ConfStat"),
  scClassify = summarize_labeling("label_scClassify_NoReject", "label_scClassify_ConfStat"),
  Seurat = summarize_labeling("label_Seurat_NoReject", "label_Seurat_ConfStat")
)

# 查看結果
print(results_list)





