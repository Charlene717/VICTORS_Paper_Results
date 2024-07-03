# Mapping and annotating query datasets
# https://satijalab.org/seurat/articles/integration_mapping

# Multimodal reference mapping
# Seurat v4 Reference Mapping
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load Packages ####
## Load packages by CRAN
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("caret")) install.packages("caret"); library(caret)

#### Load Data ####
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")

seuratObject_Sample <- pbmc.rna
seuratObject_Ref <- pbmc.rna
seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]



# 1. 準備和預處理數據
# 預處理參考數據
seuratObject_Ref <- NormalizeData(seuratObject_Ref)
seuratObject_Ref <- FindVariableFeatures(seuratObject_Ref)

# 預處理待標註數據
seuratObject_Sample <- NormalizeData(seuratObject_Sample)
seuratObject_Sample <- FindVariableFeatures(seuratObject_Sample)

# 2. 找到轉移錨點並轉移數據
# 找到轉移錨點
anchors <- FindTransferAnchors(reference = seuratObject_Ref, query = seuratObject_Sample, dims = 1:30)

# 轉移數據
predictions <- TransferData(anchorset = anchors, refdata = seuratObject_Ref$Actual_Cell_Type, dims = 1:30)

# 將預測結果添加到待標註數據中
seuratObject_Sample <- AddMetaData(seuratObject_Sample, metadata = predictions)

p1 <- DimPlot(seuratObject_Sample, group.by = "predicted.id", label = FALSE, label.size = 3) # + NoLegend()
p2 <- DimPlot(seuratObject_Sample, group.by = "seurat_annotations")
p1 + p2



# 3. Mapping QC

# # error
# # 計算映射質量得分
# seuratObject_Sample <- MappingScore(seuratObject_Sample, anchorset = anchors)

# 計算每個細胞的最大預測得分
seuratObject_Sample$mapping.score <- apply(seuratObject_Sample@meta.data[, grep("prediction.score", colnames(seuratObject_Sample@meta.data))], 1, max)



# 可視化映射質量得分
FeaturePlot(seuratObject_Sample, features = "mapping.score")

# 4. 檢查和保存結果
# 檢查標註結果
head(seuratObject_Sample@meta.data)

# 保存結果
saveRDS(seuratObject_Sample, file = "annotated_data.rds")

# 5. 評估和比較
# 假設有真實標註和預測標註
confusionMatrix(data = seuratObject_Sample$predicted.celltype, reference = seuratObject_Sample$true.celltype)









#######################################################
# Mapping and annotating query datasets
# https://satijalab.org/seurat/articles/integration_mapping

# 加載必要的包
library(Seurat)

# 假設參考數據集在 ref_data 中，查詢數據集在 query_data 中

# 參考數據集準備
ref_data <- NormalizeData(ref_data)
ref_data <- FindVariableFeatures(ref_data)
ref_data <- ScaleData(ref_data)
ref_data <- RunPCA(ref_data)

# 查詢數據集準備
query_data <- NormalizeData(query_data)
query_data <- FindVariableFeatures(query_data)
query_data <- ScaleData(query_data)
query_data <- RunPCA(query_data)

# 查詢數據集映射
anchors <- FindTransferAnchors(reference = ref_data, query = query_data, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = ref_data$celltype, dims = 1:30)
query_data <- AddMetaData(query_data, metadata = predictions)

# 結果可視化
query_data <- RunUMAP(query_data, reduction = "pca", dims = 1:30)
DimPlot(query_data, reduction = "umap", group.by = "predicted.id")

#######################################################
# Seurat v4 Reference Mapping
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

# 安裝並加載必要的包
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("SeuratDisk", quietly = TRUE)) devtools::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(SeuratDisk)

# 假設參考數據在 seuratObject_Ref 中
current_dir <- getwd()

# 刪除已存在的 .h5seurat 文件
if (file.exists(file.path(current_dir, "ref.h5seurat"))) {
  file.remove(file.path(current_dir, "ref.h5seurat"))
}

# 保存 Seurat 對象為 .h5seurat 文件
SaveH5Seurat(seuratObject_Ref, filename = file.path(current_dir, "ref.h5seurat"))

# 加載查詢數據集
# 假設查詢數據集在 seuratObject_Sample 中

# 加載參考數據集
reference <- LoadH5Seurat(file.path(current_dir, "ref.h5seurat"))

# 使用 MapQuery 將查詢數據集映射到參考數據集
query <- MapQuery(
  anchorset = reference,
  query = seuratObject_Sample,
  reference = reference,
  refdata = list(celltype = "celltype"),  # 假設 celltype 是我們要映射的註釋列
  reduction.model = "pca"
)

# 整合嵌入
query <- IntegrateEmbeddings(
  anchorset = reference,
  query = query,
  new.reduction.name = "ref.pca"
)

# 投影 UMAP
query <- ProjectUMAP(
  query = query,
  reference = reference,
  reduction.model = "umap"
)

# 可視化結果
DimPlot(query, reduction = "umap", group.by = "predicted.celltype")
