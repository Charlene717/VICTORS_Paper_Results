##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require(cowplot)) install.packages("cowplot"); library(cowplot)

source("Plot_CellAnnot_UMAP_Box.R")

#### Load dataset ####
load("D:/Dropbox/###_VUMC/##_Research/VICTORS/20231229_Figures/Export_GSE132044_MislabelB cell/20231212125506KYHDNV_Multi/20231212125506KYHDNV.RData")

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



export_folder <- "D:/Dropbox/###_VUMC/##_Research/VICTORS/20231229_Figures/Export_GSE132044_MislabelB cell/20231212125506KYHDNV_Multi/"
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



#### Fig1 ####

## Fig1B
## Create a Plot list
Fig1B_plots_UMAP <- list()
Fig1B_plots_Percent <- list()

Fig1B_plots_UMAP[["singleR_VICTORS_UMAP"]] <- plots_singleR_VICTORS_Act$UMAP_label
Fig1B_plots_Percent[["singleR_VICTORS_Percent"]] <- plots_singleR_VICTORS_Act$Percent_Stacked_Barchart
Fig1B_plots_UMAP[["scmap_VICTORS_UMAP"]] <- plots_scmap_VICTORS_Act$UMAP_label
Fig1B_plots_Percent[["scmap_VICTORS_Percent"]] <- plots_scmap_VICTORS_Act$Percent_Stacked_Barchart
Fig1B_plots_UMAP[["SCINA_VICTORS_UMAP"]] <- plots_SCINA_VICTORS_Act$UMAP_label
Fig1B_plots_Percent[["SCINA_VICTORS_Percent"]] <- plots_SCINA_VICTORS_Act$Percent_Stacked_Barchart
Fig1B_plots_UMAP[["scPred_VICTORS_UMAP"]] <- plots_scPred_VICTORS_Act$UMAP_label
Fig1B_plots_Percent[["scPred_VICTORS_Percent"]] <- plots_scPred_VICTORS_Act$Percent_Stacked_Barchart

## Combine all plots
combined_plot_UMAP <- plot_grid(plotlist = Fig1B_plots_UMAP, ncol = 4)
print(combined_plot_UMAP)

combined_plot_Percent <- plot_grid(plotlist = Fig1B_plots_Percent, ncol = 4)
print(combined_plot_Percent)

## Export PDF
pdf(paste0(export_folder, "/", export_name, "_MainResult_Fig1B1.pdf"),
    width = 16, height = 8)
print(combined_plot_UMAP)
dev.off()

pdf(paste0(export_folder, "/", export_name, "_MainResult_Fig1B2.pdf"),
    width = 20, height = 8)
print(combined_plot_Percent)

dev.off()




## Fig1A
seuratObject_Sample@meta.data$`Actual_Cell_Type` <- seuratObject_Sample@meta.data$`Actual Cell Type`
seuratObject_Sample@meta.data$`seurat_clusters` <- seuratObject_Sample@meta.data$`seurat clusters`
seuratObject_Sample@meta.data$singleR <- seuratObject_Sample@meta.data$label_singleR_DiagPara
seuratObject_Sample@meta.data$SCINA <- seuratObject_Sample@meta.data$label_SCINA_DiagPara
seuratObject_Sample@meta.data$scPred <- seuratObject_Sample@meta.data$label_scPred_DiagPara
seuratObject_Sample@meta.data$scmap <- seuratObject_Sample@meta.data$label_scmap_DiagPara


df2 <- data.frame(Actual_Cell_Type = as.character(seuratObject_Sample$Actual_Cell_Type),
                 seurat_clusters = as.character(seuratObject_Sample$seurat_clusters),
                 singleR = as.character(seuratObject_Sample@meta.data$singleR),
                 SCINA = as.character(seuratObject_Sample@meta.data$SCINA),
                 scPred = as.character(seuratObject_Sample@meta.data$scPred),
                 scmap = as.character(seuratObject_Sample@meta.data$scmap))

## singleR
plots_singleR <- Fun_Plot_UMAP_Bar(df2, seuratObject_Sample, Set_cluster = "singleR", Set_cluster_Title = "singleR",
                                   palette = "Set3", legend = FALSE, color_vector = color_Class)
plots_singleR_Act <- Fun_Plot_UMAP_Bar(df2, seuratObject_Sample, Set_cluster = "singleR", Set_cluster_Title = "singleR",
                                       Set_cluster2 = "Actual_Cell_Type",Set_cluster_Title2 = "Actual_Cell_Type",
                                       palette = "Set3", legend = FALSE, color_vector = color_Class)
plots_singleR$UMAP_label
plots_singleR_Act$UMAP_label

## scmap
plots_scmap_Act <- Fun_Plot_UMAP_Bar(df2, seuratObject_Sample, Set_cluster = "scmap", Set_cluster_Title = "scmap",
                                       Set_cluster2 = "Actual_Cell_Type",Set_cluster_Title2 = "Actual_Cell_Type",
                                       palette = "Set3", legend = FALSE, color_vector = color_Class)

## SCINA
plots_SCINA_Act <- Fun_Plot_UMAP_Bar(df2, seuratObject_Sample, Set_cluster = "SCINA", Set_cluster_Title = "SCINA",
                                     Set_cluster2 = "Actual_Cell_Type",Set_cluster_Title2 = "Actual_Cell_Type",
                                     palette = "Set3", legend = FALSE, color_vector = color_Class)

## scPred
plots_scPred_Act <- Fun_Plot_UMAP_Bar(df2, seuratObject_Sample, Set_cluster = "scPred", Set_cluster_Title = "scPred",
                                     Set_cluster2 = "Actual_Cell_Type",Set_cluster_Title2 = "Actual_Cell_Type",
                                     palette = "Set3", legend = FALSE, color_vector = color_Class)


## Create a Plot list
Fig1A_plots_UMAP <- list()
Fig1A_plots_Percent <- list()

Fig1A_plots_UMAP[["singleR_UMAP"]] <- plots_singleR_Act$UMAP_label
Fig1A_plots_Percent[["singleR_Percent"]] <- plots_singleR_Act$Percent_Stacked_Barchart
Fig1A_plots_UMAP[["scmap_UMAP"]] <- plots_scmap_Act$UMAP_label
Fig1A_plots_Percent[["scmap_Percent"]] <- plots_scmap_Act$Percent_Stacked_Barchart
Fig1A_plots_UMAP[["SCINA_UMAP"]] <- plots_SCINA_Act$UMAP_label
Fig1A_plots_Percent[["SCINA_Percent"]] <- plots_SCINA_Act$Percent_Stacked_Barchart
Fig1A_plots_UMAP[["scPred_UMAP"]] <- plots_scPred_Act$UMAP_label
Fig1A_plots_Percent[["scPred_Percent"]] <- plots_scPred_Act$Percent_Stacked_Barchart

## Combine all plots
combined_plot_UMAP_A <- plot_grid(plotlist = Fig1A_plots_UMAP, ncol = 4)
print(combined_plot_UMAP_A)

combined_plot_Percent_A <- plot_grid(plotlist = Fig1A_plots_Percent, ncol = 4)
print(combined_plot_Percent_A)

## Export PDF
pdf(paste0(export_folder, "/", export_name, "_MainResult_Fig1A1.pdf"),
    width = 16, height = 8)
print(combined_plot_UMAP_A)
dev.off()

pdf(paste0(export_folder, "/", export_name, "_MainResult_Fig1A2.pdf"),
    width = 20, height = 8)
print(combined_plot_Percent_A)

dev.off()



