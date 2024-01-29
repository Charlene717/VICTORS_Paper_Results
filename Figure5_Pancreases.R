##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Package #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require(cowplot)) install.packages("cowplot"); library(cowplot)


##### Load Data #####

## Baron
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Pancreas/IntCT_scRNAseq_Ref/BaronPancreasData.RData")
seuratObject$Sample <- seuratObject$orig.ident
p1 <- DimPlot(seuratObject, label = TRUE, repel = TRUE, group.by = "seurat_clusters") + NoLegend()
p2 <- DimPlot(seuratObject, label = TRUE, repel = TRUE, group.by = "Actual_Cell_Type") + NoLegend()
p3 <- DimPlot(seuratObject, label = FALSE, repel = TRUE, group.by = "Sample")

p1 + p2 + p3

