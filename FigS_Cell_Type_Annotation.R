#### Check Markers ####

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


##### Load Package #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require(cowplot)) install.packages("cowplot"); library(cowplot)

source("PlotFun_Beautify_UMAP_Box.R")
source("Set_plot_color.R")


#### Load Dataset ####
Dataset <- "GSE132044_PBMC_MislabelB"
load("D:/Dropbox/###_VUMC/##_Research/VICTORS/20231229_Figures/PBMC_GSE132044/Export_GSE132044_MislabelB cell/20231212125506KYHDNV_Multi/20231212125506KYHDNV.RData")

# Dataset <- "GSE132044_PBMC_MislabelNone"
# load("D:/Dropbox/###_VUMC/##_Research/VICTORS/20231229_Figures/PBMC_GSE132044/Export_GSE132044_MislabelNone/20231221065523SVLJPQ_Multi/20231221065523SVLJPQ.RData")

# Dataset <- "GSE132044_PBMC_MislabelCD16Mono"
# load("D:/Dropbox/###_VUMC/##_Research/VICTORS/20231229_Figures/PBMC_GSE132044/Export_GSE132044_MislabelCD16Mono/20231215144125WUXTQO_Multi/20231215144125WUXTQO.RData")

seuratObject <- seuratObject_Sample
seuratObject@meta.data$`Actual Cell Type` %>% unique()


#### Extract Annotation ####
seuratObject@meta.data$singleR <- seuratObject@meta.data$label_singleR
seuratObject@meta.data$scmap <- seuratObject@meta.data$label_scmap
seuratObject@meta.data$SCINA <- seuratObject@meta.data$label_SCINA
seuratObject@meta.data$scPred <- seuratObject@meta.data$label_scPred

colnames(seuratObject@meta.data) <- gsub(" ", "_",colnames(seuratObject@meta.data))

df <- data.frame(`Actual_Cell_Type` = as.character(seuratObject$`Actual_Cell_Type`),
                 singleR = as.character(seuratObject$`singleR`),
                 scmap = as.character(seuratObject$`scmap`),
                 SCINA = as.character(seuratObject$`SCINA`),
                 scPred = as.character(seuratObject$`scPred`),
                 `seurat_clusters` = as.character(seuratObject$`seurat_clusters`))


plots_Anno_CT_count1 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "singleR", Set_cluster_Title = "singleR",
                                         Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_vector)
print(plots_Anno_CT_count1$UMAP_label2 + plots_Anno_CT_count1$UMAP_label + plots_Anno_CT_count1$Grouped_Barchart + plots_Anno_CT_count1$Percent_Stacked_Barchart)

plots_Anno_CT_count2 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "scmap", Set_cluster_Title = "scmap",
                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_vector)
print(plots_Anno_CT_count2$UMAP_label2 + plots_Anno_CT_count2$UMAP_label + plots_Anno_CT_count2$Grouped_Barchart + plots_Anno_CT_count2$Percent_Stacked_Barchart)

plots_Anno_CT_count3 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "SCINA", Set_cluster_Title = "SCINA",
                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_vector)
print(plots_Anno_CT_count3$UMAP_label2 + plots_Anno_CT_count3$UMAP_label + plots_Anno_CT_count3$Grouped_Barchart + plots_Anno_CT_count3$Percent_Stacked_Barchart)

plots_Anno_CT_count4 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "scPred", Set_cluster_Title = "scPred",
                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_vector)
print(plots_Anno_CT_count4$UMAP_label2 + plots_Anno_CT_count4$UMAP_label + plots_Anno_CT_count4$Grouped_Barchart + plots_Anno_CT_count4$Percent_Stacked_Barchart)


## Combine all plots
plots_UMAP.lt <- list()
plots_UMAP.lt[["UMAP1"]] <- plots_Anno_CT_count1$UMAP_label
plots_UMAP.lt[["UMAP2"]] <- plots_Anno_CT_count2$UMAP_label
plots_UMAP.lt[["UMAP3"]] <- plots_Anno_CT_count3$UMAP_label
plots_UMAP.lt[["UMAP4"]] <- plots_Anno_CT_count4$UMAP_label

plots_Percent.lt <- list()
plots_Percent.lt[["Percent1"]] <- plots_Anno_CT_count1$Percent_Stacked_Barchart+ theme(legend.position = "none")
plots_Percent.lt[["Percent2"]] <- plots_Anno_CT_count2$Percent_Stacked_Barchart+ theme(legend.position = "none")
plots_Percent.lt[["Percent3"]] <- plots_Anno_CT_count3$Percent_Stacked_Barchart+ theme(legend.position = "none")
plots_Percent.lt[["Percent4"]] <- plots_Anno_CT_count4$Percent_Stacked_Barchart+ theme(legend.position = "none")

combined_UMAP_plot <- plot_grid(plotlist = plots_UMAP.lt, ncol = 4)
print(combined_UMAP_plot)

combined_Percent_plot <- plot_grid(plotlist = plots_Percent.lt, ncol = 4)
print(combined_Percent_plot)


## Export
Name_ExportFolder <- "Export_Annotation"

Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_Export <- paste0(Name_time_wo_micro, "_", Dataset)

if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_Annotation.pdf"),
    width = 15, height = 7)
print(combined_UMAP_plot)
print(combined_Percent_plot)
dev.off()


## Count Cell type on cluster
plots_CT_Clt_count <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "Actual_Cell_Type", Set_cluster_Title = "Actual_Cell_Type",
                                        palette = "Set3", legend = FALSE, color_vector = color_vector)
print(plots_CT_Clt_count$UMAP_cluster + plots_CT_Clt_count$UMAP_label + plots_CT_Clt_count$Grouped_Barchart + plots_CT_Clt_count$Percent_Stacked_Barchart)

plots_Anno_CT_count <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "Actual_Cell_Type", Set_cluster_Title = "Actual_Cell_Type",
                                         Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_vector)
print(plots_Anno_CT_count$UMAP_label2 + plots_Anno_CT_count$UMAP_label + plots_Anno_CT_count$Grouped_Barchart + plots_Anno_CT_count$Percent_Stacked_Barchart)


## Export
pdf(paste0(Name_ExportFolder,"/",Name_Export,"_CellTypeCount.pdf"),
    width = 12, height = 12)

print(plots_CT_Clt_count$UMAP_cluster + plots_CT_Clt_count$UMAP_label + plots_CT_Clt_count$Grouped_Barchart + plots_CT_Clt_count$Percent_Stacked_Barchart)
print(plots_Anno_CT_count$UMAP_label2 + plots_Anno_CT_count$UMAP_label + plots_Anno_CT_count$Grouped_Barchart + plots_Anno_CT_count$Percent_Stacked_Barchart)
print(plots_Anno_CT_count4$UMAP_label2 + plots_Anno_CT_count4$UMAP_label + plots_Anno_CT_count4$Grouped_Barchart + plots_Anno_CT_count4$Percent_Stacked_Barchart)

dev.off()



###################
#### singleR_Pruned ####
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("SingleR")) BiocManager::install("SingleR"); library(SingleR)

Num_Col <- seuratObject@meta.data$`Actual_Cell_Type` %>% unique() %>% length()
plots_DeltaDist <- plotDeltaDistribution(SingleR.lt, size = 2,ncol = Num_Col)
print(plots_DeltaDist)

seuratObject@meta.data$singleR_Pruned <- SingleR.lt@listData[["pruned.labels"]]
seuratObject@meta.data$singleR_Pruned[is.na(seuratObject@meta.data$singleR_Pruned)] <- "pruned"
df_singleR_Pruned <- data.frame(`Actual_Cell_Type` = as.character(seuratObject$`Actual_Cell_Type`),
                                singleR_Pruned = as.character(seuratObject$`singleR_Pruned`),
                                `seurat_clusters` = as.character(seuratObject$`seurat_clusters`))

plots_Anno_CT_count1 <- Fun_Plot_UMAP_Bar(df_singleR_Pruned, seuratObject, Set_cluster = "singleR_Pruned", Set_cluster_Title = "singleR_Pruned",
                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_vector)
print(plots_Anno_CT_count1$UMAP_label2 + plots_Anno_CT_count1$UMAP_label + plots_Anno_CT_count1$Grouped_Barchart + plots_Anno_CT_count1$Percent_Stacked_Barchart)


## Export
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_singleR_Pruned1.pdf"),
    width = 15, height = 4)
print(plots_DeltaDist)
dev.off()

pdf(paste0(Name_ExportFolder, "/", Name_Export, "_singleR_Pruned2.pdf"),
    width = 12, height = 12)
print(plots_Anno_CT_count1$UMAP_label2 + plots_Anno_CT_count1$UMAP_label + plots_Anno_CT_count1$Grouped_Barchart + plots_Anno_CT_count1$Percent_Stacked_Barchart)

dev.off()
