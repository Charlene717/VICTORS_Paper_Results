# Multimodal reference mapping

# Mapping and annotating query datasets
# https://satijalab.org/seurat/articles/integration_mapping

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
# ## Error
#
# #### Load Data ####
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Seurat_pbmcMultiome/Seurat_pbmcMultiome_Preprocessing.RData")
#
# seuratObject_Sample <- pbmc.rna
# seuratObject_Ref <- pbmc.rna
# seuratObject_Ref@meta.data[["Actual_Cell_Type"]] <- seuratObject_Ref@meta.data[["seurat_annotations"]]
#
#
# # Run Azimuth for cell type annotation
# seuratObject_Sample <- RunAzimuth(
#   query = seuratObject_Sample,
#   reference = list(
#     ref = seuratObject_Ref,
#     celltype_column = "Actual_Cell_Type"
#   ))
#
# seuratObject_Sample <- RunAzimuth(
#   query = seuratObject_Sample,
#   reference = seuratObject_Ref,
#   normalization.method = "SCT",
#   reference.assay = "refAssay",
#   query.assay = "RNA",
#   reduction = "spca",
#   dims = 1:30
# )
#
