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
# main_directory <- "D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/Export_scRNAseqPanc_20240711"

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

    # ## Add Name_FileID to meta.data
    # seuratObject_Sample@meta.data$FileID <- Name_FileID
    # save.image(file)

    ## Export MetaData
    write.table(data.frame(ID=rownames(seuratObject_Sample@meta.data), seuratObject_Sample@meta.data),
                file=paste0(main_directory,"/",Name_ExportFolder,"/",Name_Export,"_metadataSamp.tsv"),
                # file=paste0(Name_time_wo_micro,"_metadataSamp.tsv"),
                quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')

  }
}

# 顯示整合後的dataframe




#### Visualization ####





