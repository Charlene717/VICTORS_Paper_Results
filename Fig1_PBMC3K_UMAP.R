##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)

source("Plot_CellAnnot_UMAP_Box.R")

#### Load dataset ####
load("D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20231118_PBMC3K/Export_PBMC3K_MislabelB/20231119014028BQVNGA_Multi/20231119014028BQVNGA.RData")

source("FUN_Plot_Beautify_UMAP_Box.R")
source("Set_plot_color.R")

seuratObject_Sample@meta.data$`Actual_Cell_Type` <- seuratObject_Sample@meta.data$`Actual Cell Type`
seuratObject_Sample@meta.data$`seurat_clusters` <- seuratObject_Sample@meta.data$`seurat clusters`
seuratObject_Sample@meta.data$singleR_VICTORS <- seuratObject_Sample@meta.data$DiagPara_label_singleR_NoReject_SVGLRglmnet_ROC
seuratObject_Sample@meta.data$SCINA_VICTORS <- seuratObject_Sample@meta.data$DiagPara_label_SCINA_NoReject_SVGLRglmnet_ROC
seuratObject_Sample@meta.data$scPred_VICTORS <- seuratObject_Sample@meta.data$DiagPara_label_scPred_NoReject_Annot_SVGLRglmnet_ROC
seuratObject_Sample@meta.data$scmap_VICTORS <- seuratObject_Sample@meta.data$DiagPara_label_scmap_NoReject_SVGLRglmnet_ROC


seuratObject_Ref@meta.data$`Actual_Cell_Type` <- seuratObject_Ref@meta.data$`Actual Cell Type`
seuratObject_Ref@meta.data$`seurat_clusters` <- seuratObject_Ref@meta.data$`seurat clusters`

df <- data.frame(Actual_Cell_Type = as.character(seuratObject_Sample$Actual_Cell_Type),
                 seurat_clusters = as.character(seuratObject_Sample$seurat_clusters),
                 singleR_VICTORS = as.character(seuratObject_Sample@meta.data$singleR_VICTORS),
                 SCINA_VICTORS = as.character(seuratObject_Sample@meta.data$SCINA_VICTORS),
                 scPred_VICTORS = as.character(seuratObject_Sample@meta.data$scPred_VICTORS),
                 scmap_VICTORS = as.character(seuratObject_Sample@meta.data$scmap_VICTORS))

## singleR_VICTORS
plots_singleR_VICTORS <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = "singleR_VICTORS", Set_cluster_Title = "singleR_VICTORS",
                                         palette = "Set3", legend = FALSE, color_vector = color_Class)
print(plots_singleR_VICTORS$UMAP_cluster + plots_singleR_VICTORS$UMAP_label +
        plots_singleR_VICTORS$Grouped_Barchart + plots_singleR_VICTORS$Percent_Stacked_Barchart)

plots_singleR_VICTORS_Act <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = "singleR_VICTORS", Set_cluster_Title = "singleR_VICTORS",
                                               Set_cluster2 = "Actual_Cell_Type",Set_cluster_Title2 = "Actual_Cell_Type",
                                               palette = "Set3", legend = FALSE, color_vector = color_Class)
print(plots_singleR_VICTORS_Act$UMAP_cluster + plots_singleR_VICTORS_Act$UMAP_label +
        plots_singleR_VICTORS_Act$Grouped_Barchart + plots_singleR_VICTORS_Act$Percent_Stacked_Barchart)


## SCINA_VICTORS
plots_SCINA_VICTORS <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = "SCINA_VICTORS", Set_cluster_Title = "SCINA_VICTORS",
                                             palette = "Set3", legend = FALSE, color_vector = color_Class)
print(plots_SCINA_VICTORS$UMAP_cluster + plots_SCINA_VICTORS$UMAP_label +
        plots_SCINA_VICTORS$Grouped_Barchart + plots_SCINA_VICTORS$Percent_Stacked_Barchart)

plots_SCINA_VICTORS_Act <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = "SCINA_VICTORS", Set_cluster_Title = "SCINA_VICTORS",
                                               Set_cluster2 = "Actual_Cell_Type",Set_cluster_Title2 = "Actual_Cell_Type",
                                               palette = "Set3", legend = FALSE, color_vector = color_Class)
print(plots_SCINA_VICTORS_Act$UMAP_cluster + plots_SCINA_VICTORS_Act$UMAP_label +
        plots_SCINA_VICTORS_Act$Grouped_Barchart + plots_SCINA_VICTORS_Act$Percent_Stacked_Barchart)


## scPred_VICTORS
plots_scPred_VICTORS <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = "scPred_VICTORS", Set_cluster_Title = "scPred_VICTORS",
                                         palette = "Set3", legend = FALSE, color_vector = color_Class)
print(plots_scPred_VICTORS$UMAP_cluster + plots_scPred_VICTORS$UMAP_label +
        plots_scPred_VICTORS$Grouped_Barchart + plots_scPred_VICTORS$Percent_Stacked_Barchart)

plots_scPred_VICTORS_Act <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = "scPred_VICTORS", Set_cluster_Title = "scPred_VICTORS",
                                               Set_cluster2 = "Actual_Cell_Type",Set_cluster_Title2 = "Actual_Cell_Type",
                                               palette = "Set3", legend = FALSE, color_vector = color_Class)
print(plots_scPred_VICTORS_Act$UMAP_cluster + plots_scPred_VICTORS_Act$UMAP_label +
        plots_scPred_VICTORS_Act$Grouped_Barchart + plots_scPred_VICTORS_Act$Percent_Stacked_Barchart)


## scmap_VICTORS
plots_scmap_VICTORS <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = "scmap_VICTORS", Set_cluster_Title = "scmap_VICTORS",
                                         palette = "Set3", legend = FALSE, color_vector = color_Class)
print(plots_scmap_VICTORS$UMAP_cluster + plots_scmap_VICTORS$UMAP_label +
        plots_scmap_VICTORS$Grouped_Barchart + plots_scmap_VICTORS$Percent_Stacked_Barchart)

plots_scmap_VICTORS_Act <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = "scmap_VICTORS", Set_cluster_Title = "scmap_VICTORS",
                                               Set_cluster2 = "Actual_Cell_Type",Set_cluster_Title2 = "Actual_Cell_Type",
                                               palette = "Set3", legend = FALSE, color_vector = color_Class)
print(plots_scmap_VICTORS_Act$UMAP_cluster + plots_scmap_VICTORS_Act$UMAP_label +
        plots_scmap_VICTORS_Act$Grouped_Barchart + plots_scmap_VICTORS_Act$Percent_Stacked_Barchart)



export_folder <- "D:/Dropbox/##_GitHub/###_VUMC/CTAEvaluator_20231118_PBMC3K/Export_PBMC3K_MislabelB/"
export_name <- "VICTORS"
pdf(paste0(export_folder, "/", export_name, "_MainResult.pdf"),
    width = 11, height = 8)

print(plots_singleR_VICTORS$UMAP_cluster + plots_singleR_VICTORS$UMAP_label +
        plots_singleR_VICTORS$Grouped_Barchart + plots_singleR_VICTORS$Percent_Stacked_Barchart)
print(plots_singleR_VICTORS_Act$UMAP_cluster + plots_singleR_VICTORS_Act$UMAP_label +
        plots_singleR_VICTORS_Act$Grouped_Barchart + plots_singleR_VICTORS_Act$Percent_Stacked_Barchart)

print(plots_SCINA_VICTORS$UMAP_cluster + plots_SCINA_VICTORS$UMAP_label +
        plots_SCINA_VICTORS$Grouped_Barchart + plots_SCINA_VICTORS$Percent_Stacked_Barchart)
print(plots_SCINA_VICTORS_Act$UMAP_cluster + plots_SCINA_VICTORS_Act$UMAP_label +
        plots_SCINA_VICTORS_Act$Grouped_Barchart + plots_SCINA_VICTORS_Act$Percent_Stacked_Barchart)

print(plots_scPred_VICTORS$UMAP_cluster + plots_scPred_VICTORS$UMAP_label +
        plots_scPred_VICTORS$Grouped_Barchart + plots_scPred_VICTORS$Percent_Stacked_Barchart)
print(plots_scPred_VICTORS_Act$UMAP_cluster + plots_scPred_VICTORS_Act$UMAP_label +
        plots_scPred_VICTORS_Act$Grouped_Barchart + plots_scPred_VICTORS_Act$Percent_Stacked_Barchart)

print(plots_scmap_VICTORS$UMAP_cluster + plots_scmap_VICTORS$UMAP_label +
        plots_scmap_VICTORS$Grouped_Barchart + plots_scmap_VICTORS$Percent_Stacked_Barchart)
print(plots_scmap_VICTORS_Act$UMAP_cluster + plots_scmap_VICTORS_Act$UMAP_label +
        plots_scmap_VICTORS_Act$Grouped_Barchart + plots_scmap_VICTORS_Act$Percent_Stacked_Barchart)


dev.off()
