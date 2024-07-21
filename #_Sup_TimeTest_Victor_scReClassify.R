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

## Query
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Processed_Sample_OriQC/GSE132044_10xV2_Sample_CellNum1681_Seed123_Processed.RData")
seuratObject_Sample <- seuratObject ; rm(seuratObject)
ncol(seuratObject_Sample)

## Reference
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Processed_Ref_OriQC/GSE132044_10xV2A_Ref_CellNum1611_Seed123_Processed.RData")

seuratObject_Ref <- seuratObject ; rm(seuratObject)
ncol(seuratObject_Ref)
table(seuratObject_Ref$Actual_Cell_Type)

#### Sampling ####
# 定義函數來隨機抽取指定數量的細胞
sample_cells <- function(seurat_object, n_cells, seed = 123) {

  if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
  set.seed(seed)  # 設定隨機種子以確保結果可重現

  # 獲取細胞類型的計數
  cell_type_counts <- table(seurat_object$Actual_Cell_Type)

  # 計算每個細胞類型的抽取比例
  sample_proportions <- cell_type_counts / sum(cell_type_counts)

  # 計算每個細胞類型需要抽取的數目
  sample_sizes <- round(sample_proportions * n_cells)

  # 從每個細胞類型中隨機抽取對應數目的細胞
  sampled_cells <- unlist(lapply(names(sample_sizes), function(cell_type) {
    cell_indices <- WhichCells(seurat_object, ident = cell_type)
    sample(cell_indices, sample_sizes[cell_type])
  }))

  # 建立新對象包含抽取的細胞
  seurat_object_sampled <- subset(seurat_object, cells = sampled_cells)

  return(seurat_object_sampled)
}

# 使用範例
# 假設seuratObject_Ref是你的Seurat對象，並希望抽取1000個細胞
seuratObject_Ref <- sample_cells(seuratObject_Ref, 1000)

# 檢查抽取後的細胞類型分布
table(seuratObject_Ref_sampled$Actual_Cell_Type)





#### Time test ####
source("#_FUN_CellTypeAnnot.R")


## Run scReClassify
seuratObject_Sample <- run_scReClassify(seuratObject_Sample, Set_AnnoCol = "Cell_Type", Set_classifier = "svm", Set_percent = 1, Set_L = 10)


## Run VICTOR
VICTOR.lt <- VICTOR(seuratObject_Sample, seuratObject_Ref,
                    ActualCellTypeColumn = "Actual_Cell_Type",
                    AnnotCellTypeColumn = "label_singleR_NoReject",
                    seurat_version = "V5")

seuratObject_Sample <- VICTOR.lt$Query
seuratObject_Ref <- VICTOR.lt$Reference



