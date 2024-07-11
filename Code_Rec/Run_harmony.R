library(Seurat)
library(harmony)

# 假设 seuratObject_Sample 和 seuratObject_Ref 已经存在
# 合并 Seurat 对象，并为它们添加前缀
seuratObject_Sample <- RenameCells(seuratObject_Sample, add.cell.id = "Sample")
seuratObject_Ref <- RenameCells(seuratObject_Ref, add.cell.id = "Ref")

combined <- merge(seuratObject_Sample, y = seuratObject_Ref)

# 标准化数据
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)

# PCA 降维
combined <- RunPCA(combined)

# 使用 Harmony 进行批次效应校正
combined <- RunHarmony(combined, group.by.vars = "orig.ident")

# UMAP 降维和聚类
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

# 提取校正后的数据并拆分
cell_ids_sample <- colnames(combined)[grepl("^Sample_", colnames(combined))]
cell_ids_ref <- colnames(combined)[grepl("^Ref_", colnames(combined))]

# 创建新的 Seurat 对象并保留原有的 Metadata
seuratObject_Sample_harmony <- subset(combined, cells = cell_ids_sample)
seuratObject_Ref_harmony <- subset(combined, cells = cell_ids_ref)

# 创建新的 Seurat 对象并保留原有的 Metadata
seuratObject_Sample_harmony <- subset(combined, cells = cell_ids_sample)
seuratObject_Ref_harmony <- subset(combined, cells = cell_ids_ref)

# # 确认新对象的结构
# print(seuratObject_Sample_harmony)
# print(seuratObject_Ref_harmony)

DimPlot(seuratObject_Sample_harmony, reduction = "umap")
DimPlot(seuratObject_Ref_harmony, reduction = "umap")
DimPlot(combined, reduction = "umap")
DimPlot(combined, reduction = "umap",group.by = "orig.ident" )

seuratObject_Sample <- seuratObject_Sample_harmony ;rm(seuratObject_Sample_harmony)
seuratObject_Ref <- seuratObject_Ref_harmony ; rm(seuratObject_Ref_harmony)
