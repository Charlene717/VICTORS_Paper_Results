##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

if(!require('dplyr')) {install.packages('dplyr'); library(dplyr)}
if(!require('readr')) {install.packages('readr'); library(readr)}

# 讀取資料
file_path <- "D:\\Dropbox\\##_GitHub\\###_VUMC\\VICTORS_Paper_Results\\#_Revision_20240818_Fin\\Figure 3 GSE132044 Within-platform\\IntegrateAll_20240817164107OSN_GSE132044_SamePlatform.tsv"
data <- read_tsv(file_path)

# 計算VICTOR對上對應方法的Median Accuracy差異
victor_diff <- data %>%
  filter(grepl("_VICTOR", Method)) %>%  # 選擇VICTOR方法
  mutate(Original_Method = gsub("_VICTOR", "", Method)) %>%  # 去除VICTOR以取得對應的原始方法
  left_join(data, by = c("Mislabel_CellType" = "Mislabel_CellType", "Original_Method" = "Method"), suffix = c("_VICTOR", "_Original")) %>%  # 合併原始方法
  mutate(Median_Difference = Median_VICTOR - Median_Original) %>%  # 計算Median差異
  select(Mislabel_CellType, Original_Method, Median_Difference)  # 選取所需的欄位

# 顯示結果
victor_diff
