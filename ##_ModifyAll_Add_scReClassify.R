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
# main_directory <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_GSE132044_20240720/Export_GSE132044_MislabelB cell"
main_directory <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results_GSE132044_scReClassify/Export_GSE132044_20240720"
export_directory <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results_GSE132044_scReClassify/Export_GSE132044_20240720"

# 獲取所有子目錄
subdirectories <- list.dirs(main_directory, recursive = TRUE, full.names = TRUE)

# 初始化一個空的dataframe
combined_data <- data.frame()

# 遍歷所有子目錄
for (subdir in subdirectories) {
  # 獲取當前子目錄中所有_metadataSamp.tsv結尾的檔案
  RData_files <- list.files(subdir, pattern = "\\.RData$", full.names = TRUE)

  # 讀取並合併這些檔案
  for (file in RData_files) {
    load(file)

    source("#_Run_scReClassify.R")

  }
}

# 顯示整合後的dataframe
print(combined_data)




#### Visualization ####







