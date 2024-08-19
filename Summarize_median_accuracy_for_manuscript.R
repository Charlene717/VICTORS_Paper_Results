##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}
if(!require('readr')) {install.packages('readr'); library(readr)}

Set_Dataset = "GSE132044_SamePlatform"
## Load data
if(Set_Dataset == "GSE132044_SamePlatform"){
  file_path <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/#_Revision_20240818_Fin/Figure 3 GSE132044 Within-platform/"
  file_Name <- "IntegrateAll_20240817164107OSN_GSE132044_SamePlatform.tsv"
}


data <- read_tsv( paste0(file_path, file_Name))

# 計算VICTOR對上對應方法的Median Accuracy差異並轉換為百分比形式
victor_diff <- data %>%
  filter(grepl("_VICTOR", Method)) %>%  # 選擇VICTOR方法
  mutate(Original_Method = gsub("_VICTOR", "", Method)) %>%  # 去除VICTOR以取得對應的原始方法
  left_join(data, by = c("Mislabel_CellType" = "Mislabel_CellType", "Original_Method" = "Method"), suffix = c("_VICTOR", "_Original")) %>%  # 合併原始方法
  mutate(`Median_Difference(%)` = (Median_VICTOR - Median_Original) * 100) %>%  # 計算Median差異並轉為百分比
  select(Mislabel_CellType, Original_Method, `Median_Difference(%)`)  # 選取所需的欄位

# 顯示結果
victor_diff


##########

# 排除Mislabel_CellType中的None，並計算每個Original_Method中最大和最小的Median_Difference(%)
victor_diff_filtered <- victor_diff %>%
  filter(Mislabel_CellType != "None")  # 排除Mislabel_CellType為None的行

# 找出每個Original_Method中Median_Difference(%)最大的Mislabel_CellType
max_diff <- victor_diff_filtered %>%
  group_by(Original_Method) %>%
  slice(which.max(`Median_Difference(%)`)) %>%
  mutate(Difference_Type = "Max")

# 找出每個Original_Method中Median_Difference(%)最小的Mislabel_CellType
min_diff <- victor_diff_filtered %>%
  group_by(Original_Method) %>%
  slice(which.min(`Median_Difference(%)`)) %>%
  mutate(Difference_Type = "Min")

# 合併最大和最小的結果生成新的表格
final_diff <- bind_rows(max_diff, min_diff) %>%
  arrange(Original_Method, Difference_Type)

# 顯示結果
final_diff


#### Export ####
# 將 victor_diff 輸出為 TSV 檔案
write_tsv(victor_diff, paste0(file_path,"Victor_Median_Accuracy_Diff.tsv"))

# 將 final_diff 輸出為 TSV 檔案
write_tsv(final_diff,paste0(file_path,"Victor_Median_Accuracy_Diff_Range.tsv"))



