##### Visualization ####
if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)

source("PlotFun_Histogram.R")
# source("PlotFun_CombineFigs.R")
source("Set_plot_color.R")

metadata <- seuratObject_Sample@meta.data

metadata <- metadata %>%
  mutate_all(~ifelse(is.na(.), "NA", .))
metadata$seurat_clusters <- as.character(metadata$seurat_clusters)


#### Histogram ####
CombPlot_VICTORS_Hist <- function(Metadata, Set_Annotation = "Cell_Type",
                                  Set_ActualCellType = "Actual_Cell_Type",
                                  Set_Cluster = "seurat_clusters",
                                  Set_fill = "Class", # "State"
                                  Set_color = NULL) {
  if(Set_fill=="State"){
    Start_word <- "Diag_VICTORS_"
    End_word <- "_StatROC"

  }else if(Set_fill=="Class"){
    Start_word <- "DiagPara_VICTORS_"
    End_word <- "_ROC"
  }


  ## Actual Cell Type
  Plot_1 <- plot_histogram(Metadata, Set_ActualCellType, paste0(Start_word, Set_Annotation, End_word),
                           color_vector = Set_color, position_type = "stack")

  Plot_2 <- plot_histogram(Metadata, Set_ActualCellType, paste0(Start_word, Set_Annotation, End_word),
                           color_vector = Set_color, type = "proportion")

  Plot_3 <- plot_histogram(Metadata, Set_ActualCellType, paste0(Start_word, Set_Annotation, End_word),
                           color_vector = Set_color, type = "count")

  ## Cluster
  Plot_4 <- plot_histogram(Metadata, Set_Cluster, paste0(Start_word, Set_Annotation, End_word),
                           color_vector = Set_color, position_type = "stack")

  Plot_5 <- plot_histogram(Metadata,
                           Set_Cluster, paste0(Start_word, Set_Annotation, End_word),
                           color_vector = Set_color, type = "proportion")

  Plot_6 <- plot_histogram(Metadata, Set_Cluster, paste0(Start_word, Set_Annotation, End_word),
                           color_vector = Set_color, type = "count")

  ## Annotation
  Plot_7 <- plot_histogram(Metadata, Set_Annotation, paste0(Start_word, Set_Annotation, End_word),
                           color_vector = Set_color, position_type = "stack")

  Plot_8 <- plot_histogram(Metadata, Set_Annotation, paste0(Start_word, Set_Annotation, End_word),
                           color_vector = Set_color, type = "proportion")

  Plot_9 <- plot_histogram(Metadata, Set_Annotation, paste0(Start_word, Set_Annotation, End_word),
                           color_vector = Set_color, type = "count")

  ## Combine Plot
  CombPlot <- Plot_1 + Plot_2 + Plot_3 + Plot_4 + Plot_5 + Plot_6 + Plot_7 + Plot_8 + Plot_9
  return(CombPlot)
}

CombPlot <- CombPlot_VICTORS_Hist(metadata, Set_Annotation = score_types[1],
                                  Set_fill = "State", Set_color = color_State)
print(CombPlot)



#### UMAP ####

## Count CTAnnot_label on cluster
df <- data.frame(Actual_Cell_Type = as.character(seuratObject_Sample$Actual_Cell_Type),
                 seurat_clusters = as.character(seuratObject_Sample@meta.data[, "seurat_clusters"]))

for (type in score_types[1:4]) {
  df[[type]] <- as.character(seuratObject_Sample@meta.data[, type])
  df[[paste0("Diag_VICTORS_", type, "_StatROC")]] <- as.character(seuratObject_Sample@meta.data[, paste0("Diag_VICTORS_", type, "_StatROC")])
  df[[paste0("DiagPara_VICTORS_", type, "_ROC")]] <- as.character(seuratObject_Sample@meta.data[, paste0("DiagPara_VICTORS_", type, "_ROC")])
}

plots_C1 <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = paste0("DiagPara_VICTORS_", score_types[1], "_ROC"), Set_cluster_Title = paste0("DiagPara_VICTORS_", score_types[1], "_ROC"),
                            Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_Class)
plots_C2 <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = paste0("DiagPara_VICTORS_", score_types[2], "_ROC"), Set_cluster_Title = paste0("DiagPara_VICTORS_", score_types[2], "_ROC"),
                            Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_Class)
plots_C3 <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = paste0("DiagPara_VICTORS_", score_types[3], "_ROC"), Set_cluster_Title = paste0("DiagPara_VICTORS_", score_types[3], "_ROC"),
                            Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_Class)
plots_C4 <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = paste0("DiagPara_VICTORS_", score_types[4], "_ROC"), Set_cluster_Title = paste0("DiagPara_VICTORS_", score_types[4], "_ROC"),
                            Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_Class)

plot_C_UMAP <- plots_C1$UMAP_label + plots_C2$UMAP_label + plots_C3$UMAP_label + plots_C4$UMAP_label
print(plot_C_UMAP)
plot_C_Barchart <- plots_C1$Grouped_Barchart + plots_C2$Grouped_Barchart + plots_C3$Grouped_Barchart + plots_C4$Grouped_Barchart
print(plot_C_Barchart)
plot_C_Stacked_Barchart <- plots_C1$Percent_Stacked_Barchart + plots_C2$Percent_Stacked_Barchart + plots_C3$Percent_Stacked_Barchart + plots_C4$Percent_Stacked_Barchart
print(plot_C_Stacked_Barchart)


plots_1 <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = paste0("Diag_VICTORS_", score_types[1], "_StatROC"), Set_cluster_Title = paste0("Diag_VICTORS_", score_types[1], "_StatROC"),
                             Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_State)
plots_2 <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = paste0("Diag_VICTORS_", score_types[2], "_StatROC"), Set_cluster_Title = paste0("Diag_VICTORS_", score_types[2], "_StatROC"),
                             Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_State)
plots_3 <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = paste0("Diag_VICTORS_", score_types[3], "_StatROC"), Set_cluster_Title = paste0("Diag_VICTORS_", score_types[3], "_StatROC"),
                             Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_State)
plots_4 <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = paste0("Diag_VICTORS_", score_types[4], "_StatROC"), Set_cluster_Title = paste0("Diag_VICTORS_", score_types[4], "_StatROC"),
                             Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_State)

plot_UMAP <- plots_1$UMAP_label + plots_2$UMAP_label + plots_3$UMAP_label + plots_4$UMAP_label
print(plot_UMAP)
plot_Barchart <- plots_1$Grouped_Barchart + plots_2$Grouped_Barchart + plots_3$Grouped_Barchart + plots_4$Grouped_Barchart
print(plot_Barchart)
plot_Stacked_Barchart <- plots_1$Percent_Stacked_Barchart + plots_2$Percent_Stacked_Barchart + plots_3$Percent_Stacked_Barchart + plots_4$Percent_Stacked_Barchart
print(plot_Stacked_Barchart)



# print(plot_1$UMAP_label + plot_2$UMAP_label + plot_3$UMAP_label + plot_4$UMAP_label)
# # print(plot_1$UMAP_label2 + plot_1$UMAP_label + plot_1$Grouped_Barchart + plot_1$Percent_Stacked_Barchart)


#### Box Plot ####





#### Export ####
## PDF
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_AnnoDiagnosis_Hist.pdf"), width = 17, height = 17)
for (score_type in score_types) {
  CombPlot_Class <- CombPlot_VICTORS_Hist(metadata, Set_Annotation = score_type,
                                          Set_fill = "Class", Set_color = color_Class)
  print(CombPlot_Class)

  CombPlot_State <- CombPlot_VICTORS_Hist(metadata, Set_Annotation = score_type,
                                          Set_fill = "State", Set_color = color_State)
  print(CombPlot_State)

}

dev.off()


pdf(paste0(Name_ExportFolder, "/", Name_Export, "_AnnoDiagnosis_Main.pdf"), width = 12, height = 12)
print(plot_C_UMAP)
print(plot_C_Barchart)
# print(plot_C_Stacked_Barchart)
print(plot_UMAP)
print(plot_Barchart)
# print(plot_Stacked_Barchart)

dev.off()
