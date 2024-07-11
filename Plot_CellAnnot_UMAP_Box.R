plot_seurat_data <- function(seurat_obj_Sample, seurat_obj_Ref ,
                             label_column_name, label_column_name2 = NULL,
                             DiagPara_column_name = NULL,
                             cluster_column_name = "seurat_clusters",
                             export_folder = "", export_name =""
                             ) {

  # source("Load_Marker_Gene.R")
  # Marker_gene.vector <- unlist(Marker_gene_pbmc.list)
  # # FeaturePlot(seurat_obj_Sample, features = Marker_gene.vector)

  source("FUN_Plot_Beautify_UMAP_Box.R")
  source("Set_plot_color.R")

  df <- data.frame(Actual_Cell_Type = as.character(seurat_obj_Sample$Actual_Cell_Type))
  df[[label_column_name]] <- as.character(seurat_obj_Sample@meta.data[, label_column_name])
  df[[cluster_column_name]] <- as.character(seurat_obj_Sample@meta.data[, cluster_column_name])

  df_Ref <- data.frame(Actual_Cell_Type = as.character(seurat_obj_Ref$Actual_Cell_Type))
  df_Ref[[cluster_column_name]] <- as.character(seurat_obj_Ref@meta.data[, cluster_column_name])


  ## Count CTAnnot_label on cluster
  plots_CTAnnot_Clt_count <- Fun_Plot_UMAP_Bar(df, seurat_obj_Sample, Set_cluster = label_column_name, Set_cluster_Title = label_column_name,
                                               palette = "Set3", legend = FALSE, color_vector = color_CellType)
  print(plots_CTAnnot_Clt_count$UMAP_cluster + plots_CTAnnot_Clt_count$UMAP_label + plots_CTAnnot_Clt_count$Grouped_Barchart + plots_CTAnnot_Clt_count$Percent_Stacked_Barchart)


  # Count CTAnnot_label on Cell Type
  plots_CTAnnot_CT_count <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = label_column_name, Set_cluster_Title = label_column_name,
                                              Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_CellType)
  print(plots_CTAnnot_CT_count$UMAP_label2 + plots_CTAnnot_CT_count$UMAP_label + plots_CTAnnot_CT_count$Grouped_Barchart + plots_CTAnnot_CT_count$Percent_Stacked_Barchart)



  if (!is.null(label_column_name2) && label_column_name2 %in% colnames(seurat_obj_Sample@meta.data)) {
    df[[label_column_name2]] <- as.character(seurat_obj_Sample@meta.data[, label_column_name2])

    ## Count CTAnnot_label on cluster
    plots_CTAnnot_Clt_count2 <- Fun_Plot_UMAP_Bar(df, seurat_obj_Sample, Set_cluster = label_column_name2, Set_cluster_Title = label_column_name2,
                                                 palette = "Set3", legend = FALSE, color_vector = color_CellType)
    print(plots_CTAnnot_Clt_count2$UMAP_cluster + plots_CTAnnot_Clt_count2$UMAP_label + plots_CTAnnot_Clt_count2$Grouped_Barchart + plots_CTAnnot_Clt_count2$Percent_Stacked_Barchart)

    # Count CTAnnot_label on Cell Type
    plots_CTAnnot_CT_count2 <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = label_column_name2, Set_cluster_Title = label_column_name2,
                                                Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_CellType)
    print(plots_CTAnnot_CT_count2$UMAP_label2 + plots_CTAnnot_CT_count2$UMAP_label + plots_CTAnnot_CT_count2$Grouped_Barchart + plots_CTAnnot_CT_count2$Percent_Stacked_Barchart)
  }


  if (!is.null(DiagPara_column_name) && DiagPara_column_name %in% colnames(seurat_obj_Sample@meta.data)) {
    df[[DiagPara_column_name]] <- as.character(seurat_obj_Sample@meta.data[, DiagPara_column_name])

    ## Count CTAnnot_label on cluster
    plots_CTAnnot_Clt_count_Diag <- Fun_Plot_UMAP_Bar(df, seurat_obj_Sample, Set_cluster = DiagPara_column_name, Set_cluster_Title = DiagPara_column_name,
                                                  palette = "Set3", legend = FALSE, color_vector = color_Class)
    print(plots_CTAnnot_Clt_count2$UMAP_cluster + plots_CTAnnot_Clt_count_Diag$UMAP_label + plots_CTAnnot_Clt_count_Diag$Grouped_Barchart + plots_CTAnnot_Clt_count_Diag$Percent_Stacked_Barchart)

    # Count CTAnnot_label on Cell Type
    plots_CTAnnot_CT_count_Diag <- Fun_Plot_UMAP_Bar(df, seuratObject_Sample, Set_cluster = DiagPara_column_name, Set_cluster_Title = DiagPara_column_name,
                                                 Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_Class)
    print(plots_CTAnnot_CT_count2$UMAP_label2 + plots_CTAnnot_CT_count_Diag$UMAP_label + plots_CTAnnot_CT_count_Diag$Grouped_Barchart + plots_CTAnnot_CT_count_Diag$Percent_Stacked_Barchart)
  }


  ## Plot count
  # Count Cell type on cluster
  plots_CT_Clt_count <- Fun_Plot_UMAP_Bar(df, seurat_obj_Sample, Set_cluster = "Actual_Cell_Type", Set_cluster_Title = "Actual_Cell_Type",
                                          palette = "Set3", legend = FALSE, color_vector = color_CellType)
  print(plots_CT_Clt_count$UMAP_cluster + plots_CT_Clt_count$UMAP_label + plots_CT_Clt_count$Grouped_Barchart + plots_CT_Clt_count$Percent_Stacked_Barchart)

  plots_CT_Clt_count_Ref <- Fun_Plot_UMAP_Bar(df_Ref, seurat_obj_Ref, Set_cluster = "Actual_Cell_Type", Set_cluster_Title = "Actual_Cell_Type(Ref)",
                                              palette = "Set3", legend = FALSE, color_vector = color_CellType)
  print(plots_CT_Clt_count_Ref$UMAP_cluster + plots_CT_Clt_count_Ref$UMAP_label + plots_CT_Clt_count_Ref$Grouped_Barchart + plots_CT_Clt_count_Ref$Percent_Stacked_Barchart)




  pdf(paste0(export_folder, "/", export_name, "_MainResult.pdf"),
      width = 11, height = 8)
  if (!is.null(DiagPara_column_name) && DiagPara_column_name %in% colnames(seurat_obj_Sample@meta.data)) {
    print(plots_CTAnnot_CT_count2$UMAP_label2 + plots_CTAnnot_CT_count_Diag$UMAP_label + plots_CTAnnot_CT_count_Diag$Grouped_Barchart + plots_CTAnnot_CT_count_Diag$Percent_Stacked_Barchart)
    print(plots_CTAnnot_Clt_count_Diag$UMAP_cluster + plots_CTAnnot_Clt_count_Diag$UMAP_label + plots_CTAnnot_Clt_count_Diag$Grouped_Barchart + plots_CTAnnot_Clt_count_Diag$Percent_Stacked_Barchart)
  }

  print(plots_CTAnnot_CT_count$UMAP_label2 + plots_CTAnnot_CT_count$UMAP_label + plots_CTAnnot_CT_count$Grouped_Barchart + plots_CTAnnot_CT_count$Percent_Stacked_Barchart)
  if (!is.null(label_column_name2) && label_column_name2 %in% colnames(seurat_obj_Sample@meta.data)) {
    print(plots_CTAnnot_CT_count2$UMAP_label2 + plots_CTAnnot_CT_count2$UMAP_label + plots_CTAnnot_CT_count2$Grouped_Barchart + plots_CTAnnot_CT_count2$Percent_Stacked_Barchart)
  }

  print(plots_CTAnnot_Clt_count$UMAP_cluster + plots_CTAnnot_Clt_count$UMAP_label + plots_CTAnnot_Clt_count$Grouped_Barchart + plots_CTAnnot_Clt_count$Percent_Stacked_Barchart)
  if (!is.null(label_column_name2) && label_column_name2 %in% colnames(seurat_obj_Sample@meta.data)) {
    print(plots_CTAnnot_Clt_count2$UMAP_cluster + plots_CTAnnot_Clt_count2$UMAP_label + plots_CTAnnot_Clt_count2$Grouped_Barchart + plots_CTAnnot_Clt_count2$Percent_Stacked_Barchart)
  }

  print(plots_CT_Clt_count$UMAP_cluster + plots_CT_Clt_count$UMAP_label + plots_CT_Clt_count$Grouped_Barchart + plots_CT_Clt_count$Percent_Stacked_Barchart)
  print(plots_CT_Clt_count_Ref$UMAP_cluster + plots_CT_Clt_count_Ref$UMAP_label + plots_CT_Clt_count_Ref$Grouped_Barchart + plots_CT_Clt_count_Ref$Percent_Stacked_Barchart)

  # try({ print(FeaturePlot(seurat_obj_Sample, features = Marker_gene.vector)) })

  dev.off()
}

# ## Test Function
# plot_seurat_data(seuratObject_Sample, seuratObject_Ref,
#                  label_column_name = "label_singleR",label_column_name = "label_singleR_NoRejection",
#                  cluster_column_name = "seurat_clusters",
#                  export_folder = Name_ExportFolder,
#                  export_name = Name_Export)
