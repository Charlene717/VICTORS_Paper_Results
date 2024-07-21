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
    cell_indices <- which(seurat_object$Actual_Cell_Type == cell_type)
    sample(cell_indices, sample_sizes[cell_type])
  }))

  # 建立新對象包含抽取的細胞
  seurat_object_sampled <- subset(seurat_object, cells = sampled_cells)

  return(seurat_object_sampled)
}


# 使用範例
# 假設seuratObject_Ref是你的Seurat對象，並希望抽取1000個細胞
seuratObject_Ref_sampled <- sample_cells(seuratObject_Ref, 1000)

# 檢查抽取後的細胞類型分布
table(seuratObject_Ref_sampled$Actual_Cell_Type)

#### Time test ####
source("#_FUN_CellTypeAnnot.R")

# 設定不同細胞數目
cell_numbers <- c(200, 400, 600, 800, 1000, 1200, 1400, 1600)

# 初始化時間結果
time_results <- data.frame(Cell_Number = cell_numbers, scReClassify_Time = NA, VICTOR_Time = NA)

# 測試不同seuratObject_Sample數目運行scReClassify和VICTOR需要多少時間
for (i in 1:length(cell_numbers)) {
  n_cells <- cell_numbers[i]
  seuratObject_Sample_subset <- sample_cells(seuratObject_Sample, n_cells)

  # 測試scReClassify時間
  start_time <- Sys.time()
  seuratObject_Sample_subset <- Fun_scReClassify(seuratObject_Sample_subset, Set_AnnoCol = "Cell_Type", Set_classifier = "svm", Set_percent = 1, Set_L = 10)
  end_time <- Sys.time()
  time_results$scReClassify_Time[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # 測試VICTOR時間
  start_time <- Sys.time()
  VICTOR.lt <- VICTOR(seuratObject_Sample_subset, seuratObject_Ref_sampled,
                      ActualCellTypeColumn = "Actual_Cell_Type",
                      AnnotCellTypeColumn = "Cell_Type",
                      seurat_version = "V5")
  end_time <- Sys.time()
  time_results$VICTOR_Time[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))

  seuratObject_Sample_subset <- VICTOR.lt$Query
  # seuratObject_Ref_sampled <- VICTOR.lt$Reference
}

#### 視覺化 ####
library(ggplot2)

ggplot(time_results, aes(x = Cell_Number)) +
  geom_line(aes(y = scReClassify_Time, color = "scReClassify")) +
  geom_point(aes(y = scReClassify_Time, color = "scReClassify")) +
  geom_text(aes(y = scReClassify_Time, label = round(scReClassify_Time, 1), color = "scReClassify"), vjust = -0.5, size = 5) +
  geom_line(aes(y = VICTOR_Time, color = "VICTOR")) +
  geom_point(aes(y = VICTOR_Time, color = "VICTOR")) +
  geom_text(aes(y = VICTOR_Time, label = round(VICTOR_Time, 1), color = "VICTOR"), vjust = -0.5, size = 5) +
  labs(title = "Time taken by scReClassify and VICTOR for different sample sizes",
       x = "Number of Cells in seuratObject_Sample",
       y = "Time (seconds)",
       color = "Method") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15)
  ) +
  scale_color_manual(values = c("scReClassify" = "#7373B9", "VICTOR" = "#FF60AF"))
