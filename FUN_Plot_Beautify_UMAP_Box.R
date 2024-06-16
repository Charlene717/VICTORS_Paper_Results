## Load Packages
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)

FUN_Plot_Beautify_UMAP <- function(UMAP){
  UMAP <- UMAP+theme(panel.background = element_rect(fill = NA, color = "black", linewidth = 1.2)) +
    coord_fixed()

  return(UMAP)
}


# ## Test the function
# ## Check original cell type
# plt.UMAP_CTOri <- DimPlot(seuratObject, reduction = "umap", group.by ="Cell_Type" ,label = FALSE, pt.size = 1)# + NoLegend()
# plt.UMAP_CTOri %>% FUN_Plot_Beautify_UMAP()

###################################################################################################################
# Function to generate color_vector
Fun_Generate_Color_Vector <- function(unique_items, palette = "Set3", provided_color_vector = NULL){
  library(RColorBrewer)

  # If the user has provided a Fin_color_vector, use that, else generate Fin_color_vector using palette
  if (!is.null(provided_color_vector)) {
    # get missing colors and fill them using colorRampPalette or brewer.pal
    missing_colors <- setdiff(names(provided_color_vector), unique_items)
    if(length(missing_colors) > 0) {
      if(length(missing_colors) > brewer.pal.info[palette, "maxcolors"]) {
        Fin_color_vector_missing <- colorRampPalette(brewer.pal(brewer.pal.info[palette, "maxcolors"], palette))(length(missing_colors))
      } else {
        Fin_color_vector_missing <- brewer.pal(length(missing_colors), palette)
      }
      Fin_color_vector <- c(provided_color_vector, setNames(Fin_color_vector_missing, missing_colors))
    }
  } else {
    if(length(unique_items) > brewer.pal.info[palette, "maxcolors"]) {
      Fin_color_vector <- colorRampPalette(brewer.pal(brewer.pal.info[palette, "maxcolors"], palette))(length(unique_items))
    } else {
      Fin_color_vector <- brewer.pal(length(unique_items), palette)
    }
    Fin_color_vector <- setNames(Fin_color_vector, unique_items)
  }

  return(Fin_color_vector)
}


###################################################################################################################
## Load Packages
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("RColorBrewer")) install.packages("RColorBrewer"); library(RColorBrewer)

Fun_Plot_UMAP_Box <- function(df, seuratObject, Set_cluster = "Cell_Type", Set_cluster_Title = "Cell_Type",
                              Set_score = "SingleR.scores", xlab = "", ylab = "Score",
                              palette = "Set3", legend = TRUE, color_vector = NULL) {
  library(RColorBrewer)
  library(ggplot2)
  library(dplyr)

  # Get the number of unique cell types
  unique_celltypes <- unique(df[[Set_cluster]])

  # Generate color_vector
  color_mapping <- Fun_Generate_Color_Vector(unique_celltypes, palette, color_vector)

  # # Create a mapping from cell type to color
  # color_mapping <- color_vector

  # Merge the color mapping into the original dataframe
  df$Color <- color_mapping[df[[Set_cluster]]]

  # Generate boxplot with colors from UMAP
  plt_Box <- ggplot(df, aes_string(x = Set_cluster, y = Set_score, fill = "Color")) +
    geom_boxplot() +
    theme_light() +
    labs(title = Set_cluster_Title)+
    theme(plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
          plot.title.position = "panel",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Rotate x-axis text and increase size
          axis.text.y = element_text(size = 12),   # Increase the size of the y-axis tick labels
          axis.title.x = element_text(size = 14),  # Increase the size of the x-axis title
          axis.title.y = element_text(size = 14),  # Increase the size of the y-axis title
          panel.grid.major = element_line(color = "grey", linetype = "dashed"),
          panel.grid.minor = element_line(color = "grey", linetype = "dashed")) +
    theme(aspect.ratio = 1, panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2)) +
    labs(x = xlab, y = ylab) +
    scale_fill_identity() +
    coord_cartesian(ylim = c(0, 1))


  # Generate UMAP plot with same colors
  umap_df <- as.data.frame(Embeddings(seuratObject[["umap"]]))
  umap_df[[Set_cluster]] <- unlist(seuratObject[[Set_cluster]])
  umap_df$Color <- color_mapping[umap_df[[Set_cluster]]]

  plt_UMAP_label <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = !!sym(Set_cluster), fill = !!sym(Set_cluster)), alpha = 0.7) +
    theme_classic() +
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    theme(panel.background = element_rect(fill = NA, color = "black", linewidth = 1.2)) +
    coord_fixed() +
    labs(title = Set_cluster_Title,
         x = "UMAP 1",  # adjust as needed
         y = "UMAP 2") + # adjust as needed
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          plot.title.position = "panel",
          axis.title.x = element_text(size = 14),  # Increase the size of the x-axis title
          axis.title.y = element_text(size = 14),  # Increase the size of the y-axis title
          axis.text.x = element_text(size = 12),   # Increase the size of the x-axis tick labels
          axis.text.y = element_text(size = 12),  # Increase the size of the y-axis tick labels
          aspect.ratio = 1)

  # If legend is requested, place it on the right
  if (legend) {
    plt_UMAP_label <- plt_UMAP_label + theme(legend.position = "right")
  } else {
    plt_UMAP_label <- plt_UMAP_label + theme(legend.position = "none")
  }

  plt_UMAP_Score <- FeaturePlot(seuratObject, features = Set_score)+
    theme(panel.background = element_rect(fill = NA, color = "black", linewidth = 1.2)) +
    coord_fixed()

  return(list("Boxplot" = plt_Box, "UMAP_label" = plt_UMAP_label, "UMAP_Score" = plt_UMAP_Score))
}



# ## Test the function
# ## Add SingleR labels to Seurat object metadata
# seuratObject@meta.data[["singleR_label"]] <- SingleR.lt$labels
# # seuratObject@meta.data[["singleR_label"]] <- SingleR.lt$first.labels
# ## seuratObject@meta.data[["singleR_label"]] <- SingleR.lt$pruned.labels
#
# # # Plot singleR_label
# # plt.UMAP_CTPred <- DimPlot(seuratObject, reduction = "umap", group.by = paste0("singleR_label") ,label = TRUE, pt.size = 0.5)# + NoLegend()
# # plt.UMAP_CTPred
#
#
# ## Add SingleR scores to Seurat object metadata
# labels_SingleR <- SingleR.lt$labels  # labels_SingleR <- seuratObject@meta.data[, paste0("singleR_", Name_Note)]
# seuratObject[["SingleR.scores"]] <- sapply(1:length(labels_SingleR), function(i) SingleR.lt$scores[i, labels_SingleR[i]])
# # seuratObject[["SingleR.tuning.scores"]] <- SingleR.lt@listData[["tuning.scores"]]@listData[["first"]]
#
# #Bug# seuratObject[["SingleR.scores"]] <- SingleR.lt$scores[,1]  # assuming scores is the first column
#
# # # Plot SingleR scores
# # FeaturePlot(seuratObject, features = "SingleR.scores")
#
# # Store the integrated metric into the metadata of the Seurat object
# seuratObject@meta.data[["integrated_metric"]] <- integrated_metric
# FeaturePlot(seuratObject, features = "integrated_metric")
#
#
# ## Set color
# color_vector <- list(
#   "Memory CD4 T" = "#1f77b4",  # blue
#   "B" = "#ff7f0e",             # orange
#   "CD14+ Mono" = "#2ca02c",    # green
#   "NK" = "#d62728",            # red
#   "CD8 T" = "#9467bd",         # purple
#   "Naive CD4 T" = "#8c564b",   # brown
#   "FCGR3A+ Mono" = "#e377c2",  # pink
#   "Unknown" = "#7f7f7f",       # grey
#   "DC" = "#bcbd22",            # olive
#   "Platelet" = "#17becf"       # light blue
# )
#
# df <- data.frame(Cell_Type = as.character(seuratObject$Cell_Type),
#                  singleR_label = as.character(seuratObject$singleR_label),
#                  seurat_clusters = as.character(seuratObject$seurat_clusters),
#                  SingleR.scores = seuratObject$SingleR.scores,
#                  integrated_metric = seuratObject$integrated_metric)
#
# ## SingleR.scores & singleR_label
# plots_singleR_label <- Fun_Plot_UMAP_Box(df, seuratObject, Set_cluster = "singleR_label",
#                            Set_score = "SingleR.scores", ylab = "SingleR.scores",
#                            color_vector = color_vector,legend = FALSE)
# print(plots_singleR_label$UMAP_label + plots_singleR_label$Boxplot + plots_singleR_label$UMAP_Score)
#
# ## SingleR.scores & seurat_clusters
# plots_seurat_clusters <- Fun_Plot_UMAP_Box(df, seuratObject, Set_cluster = "seurat_clusters",
#                                          Set_score = "SingleR.scores", ylab = "SingleR.scores",
#                                          legend = FALSE)
# print(plots_seurat_clusters$UMAP_label + plots_seurat_clusters$Boxplot + plots_seurat_clusters$UMAP_Score)
#
# ## SingleR.scores & Cell_Type
# plots_Ori_label <- Fun_Plot_UMAP_Box(df, seuratObject, Set_cluster = "Cell_Type",
#                                          Set_score = "SingleR.scores", ylab = "SingleR.scores",
#                                          color_vector = color_vector,legend = FALSE)
# print(plots_Ori_label$UMAP_label + plots_Ori_label$Boxplot + plots_Ori_label$UMAP_Score)
# FeaturePlot(seuratObject, features = "PPBP") %>% FUN_Plot_Beautify_UMAP()
#
#
# ## integrated_metric & singleR_label
# plots_singleR_label <- Fun_Plot_UMAP_Box(df, seuratObject, Set_cluster = "singleR_label",
#                                          Set_score = "integrated_metric", ylab = "integrated_metric",
#                                          color_vector = color_vector,legend = FALSE)
# print(plots_singleR_label$UMAP_label + plots_singleR_label$Boxplot + plots_singleR_label$UMAP_Score)
#
# ## integrated_metric & seurat_clusters
# plots_seurat_clusters <- Fun_Plot_UMAP_Box(df, seuratObject, Set_cluster = "seurat_clusters",
#                                            Set_score = "integrated_metric", ylab = "integrated_metric",
#                                            legend = FALSE)
# print(plots_seurat_clusters$UMAP_label + plots_seurat_clusters$Boxplot + plots_seurat_clusters$UMAP_Score)
#


####################################################################################
Fun_Plot_UMAP_Bar <- function(df, seuratObject, Set_cluster = "Cell_Type", Set_cluster_Title = "Cell_Type",
                              Set_cluster2 = "seurat_clusters",Set_cluster_Title2 = "seurat_clusters",
                              palette = "Set3", legend = TRUE, color_vector = NULL) { # color_vector2 = NULL
  library(ggplot2)
  library(dplyr)

  # Get the number of unique cell types
  unique_celltypes <- unique(df[[Set_cluster]])

  # Generate color_vector
  color_mapping <- Fun_Generate_Color_Vector(unique_celltypes, palette, color_vector)

  # # Create a mapping from cell type to color
  # color_mapping <- color_vector

  # Merge the color mapping into the original dataframe
  df$Color <- color_mapping[df[[Set_cluster]]]

  # Generate UMAP plot with same colors
  umap_df <- as.data.frame(Embeddings(seuratObject[["umap"]]))
  umap_df[[Set_cluster]] <- unlist(seuratObject[[Set_cluster]])
  umap_df[[Set_cluster2]] <- unlist(seuratObject[[Set_cluster2]])
  umap_df[["seurat_clusters"]] <- unlist(seuratObject[["seurat_clusters"]])
  umap_df$Color <- color_mapping[umap_df[[Set_cluster]]]

  plt_UMAP_label <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = !!sym(Set_cluster), fill = !!sym(Set_cluster)), alpha = 0.7) +
    theme_classic() +
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    theme(panel.background = element_rect(fill = NA, color = "black", linewidth = 1.2)) +
    coord_fixed() +
    labs(title = paste0(Set_cluster_Title),
         x = "UMAP 1",
         y = "UMAP 2") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          plot.title.position = "panel",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          aspect.ratio = 1)

  plt_UMAP_label2 <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = !!sym(Set_cluster2), fill = !!sym(Set_cluster2)), alpha = 0.7) +
    theme_classic() +
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    theme(panel.background = element_rect(fill = NA, color = "black", linewidth = 1.2)) +
    coord_fixed() +
    labs(title = paste0(Set_cluster_Title2),
         x = "UMAP 1",
         y = "UMAP 2") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          plot.title.position = "panel",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          aspect.ratio = 1)

  plt_UMAP_cluster <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = !!sym("seurat_clusters"), fill = !!sym("seurat_clusters")), alpha = 0.7) +
    theme_classic() +
    scale_color_manual(values = Fun_Generate_Color_Vector(unique(df[["seurat_clusters"]]), palette = "Set3", provided_color_vector = NULL)) +
    scale_fill_manual(values = Fun_Generate_Color_Vector(unique(df[["seurat_clusters"]]), palette = "Set3", provided_color_vector = NULL)) +
    theme(panel.background = element_rect(fill = NA, color = "black", linewidth = 1.2)) +
    coord_fixed() +
    labs(title = paste0("seurat_clusters"),
         x = "UMAP 1",
         y = "UMAP 2") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
          plot.title.position = "panel",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          aspect.ratio = 1)

  if (!legend) {
    plt_UMAP_label <- plt_UMAP_label + theme(legend.position = "none")
    plt_UMAP_label2 <- plt_UMAP_label2 + theme(legend.position = "none")
  }

  # Generate the Grouped barchart
  group_barchart_df <- as.data.frame(table(seuratObject@meta.data[[Set_cluster2]],
                                           seuratObject@meta.data[[Set_cluster]]))

  plt_Grouped_Barchart <- ggplot(group_barchart_df, aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = color_vector) +
    theme_minimal() +
    labs(title = paste0(Set_cluster_Title," on ",Set_cluster2))+
    theme(plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
          plot.title.position = "panel",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_vline(data = group_barchart_df, aes(xintercept = as.numeric(Var1) + 0.5),
               linetype = "dashed", color = "grey") +
    theme(panel.background = element_rect(fill = NA, color = "black", linewidth = 1.2)) +
    labs(x = Set_cluster2, y = "Cell count", fill = Set_cluster) +
    theme(plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
          plot.title.position = "panel",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          aspect.ratio = 1)

  if (!legend) {
    plt_Grouped_Barchart <- plt_Grouped_Barchart + theme(legend.position = "none")
  }

  # Generate the Percent stacked barchart
  group_barchart_df$Percentage <- group_barchart_df$Freq / sum(group_barchart_df$Freq)

  plt_Percent_Stacked_Barchart <- ggplot(group_barchart_df, aes(x = Var1, y = Percentage, fill = Var2)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = color_vector) +
    theme_minimal() +
    labs(title = paste0(Set_cluster_Title," on ",Set_cluster2))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = Set_cluster2, y = "Cell Percentage", fill = Set_cluster) +
    theme(panel.background = element_rect(fill = NA, color = "black", linewidth = 1.2)) +
    theme(plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
          plot.title.position = "panel",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          aspect.ratio = 1)

  return(list("UMAP_label" = plt_UMAP_label, "UMAP_label2" = plt_UMAP_label2,
              "UMAP_cluster" = plt_UMAP_cluster,
              "Grouped_Barchart" = plt_Grouped_Barchart,
              "Percent_Stacked_Barchart" = plt_Percent_Stacked_Barchart))
}



# ## Test Function
# plots_CT_count <- Fun_Plot_UMAP_Bar(df = singleR.df, seuratObject, Set_cluster = "Cell_Type", Set_cluster_Title = "Cell_Type",
#                          palette = "Set3", legend = TRUE, color_vector = NULL)
# plots_CT_count <- Fun_Plot_UMAP_Bar(df = singleR.df, seuratObject, Set_cluster = "Cell_Type", Set_cluster_Title = "Cell_Type",
#                          palette = "Set3", legend = TRUE, color_vector = color_vector)
# plots_CT_count$UMAP_label + plots_CT_count$Grouped_Barchart + plots_CT_count$Percent_Stacked_Barchart
#
#
#
# plots_CT_count <- Fun_Plot_UMAP_Bar(df = singleR.df, seuratObject, Set_cluster = "singleR_label", Set_cluster_Title = "Cell_Type",
#                          palette = "Set3", legend = FALSE, color_vector = color_vector)
# plots_CT_count$UMAP_label + plots_CT_count$Grouped_Barchart + plots_CT_count$Percent_Stacked_Barchart
