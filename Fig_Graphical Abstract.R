##### Presetting ######
rm(list = ls()) # Clean variable ##* Comment out if Run All
memory.limit(150000)


#### Load package ####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)


#### Load data ####
load("D:/Dropbox/###_VUMC/##_Research/VICTORS/20231229_Figures/PBMC_GSE132044/Export_GSE132044_MislabelB cell/20231212125506KYHDNV_Multi/20231212125506KYHDNV.RData")


## seuratObject_Sample
seuratObject_Sample@meta.data$VICTORS <- seuratObject_Sample@meta.data$Diag_SVGLRglmnet_StatROC
# Directly modifying the dataframe in your Seurat object
seuratObject_Sample@meta.data$VICTORS <- ifelse(seuratObject_Sample@meta.data$VICTORS == "T", "Reliable", "Unreliable")

seuratObject_Sample@meta.data$Query <- "Unknown"
seuratObject_Sample@meta.data$Annotation <- seuratObject_Sample@meta.data$label_singleR_NoReject
seuratObject_Sample@meta.data$`Actual_Cell_Type` <- seuratObject_Sample@meta.data$`Actual Cell Type`


# Define colors with transparency
colors_with_transparency <- c("Reliable" = "#3b4da8AA", "Unreliable" = "#bd1539AA") # AA at the end for transparency
colors_with_transparency <- c("Reliable" = "#3b4da855", "Unreliable" = "#bd153955",
                              "Unknown" = "#5b5b5c55") # 55 at the end for higher transparency

Plot_UMAP_VICTORS <- DimPlot(seuratObject_Sample, reduction = "umap", group.by = "VICTORS", label = FALSE, pt.size = 3) +
  scale_color_manual(values = colors_with_transparency) +
  NoLegend() +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
        text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        aspect.ratio = 1)

Plot_UMAP_VICTORS


Plot_UMAP_Unknown <- DimPlot(seuratObject_Sample, reduction = "umap", group.by = "Query", label = FALSE, pt.size = 3) +
  scale_color_manual(values = colors_with_transparency) +
  NoLegend() +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
        text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        aspect.ratio = 1)

Plot_UMAP_Unknown


color_CellType <- c("CD4+ T cell" = "#1f77b455",  # blue
                    "B cell" = "#ff7f0e55",             # orange
                    "CD14+ monocyte" = "#2ca02c55",    # green
                    "Natural killer cell" = "#d6272855",            # red
                    "Cytotoxic T cell" = "#9467bd55",         # purple
                    "Megakaryocyte" = "#8c564b55",   # brown
                    "CD16+ monocyte" = "#e377c255",  # pink
                    "Unknown" = "#7f7f7f55",       # grey
                    "unknown" = "#7f7f7f55",       # grey
                    "None" ="#7f7f7f55",
                    "Unassigned" = "#84e8d455",
                    "Unassign" = "#84e8d455",
                    "pruned" = "#f0ad8d55",  #"#4a1919",       # dark grey
                    "Dendritic cell" = "#bcbd2255",            # olive
                    "Plasmacytoid dendritic cell" = "#17becf55"      # light blue
)


# Plot_UMAP_Annotation <- DimPlot(seuratObject_Sample, reduction = "umap", group.by = "Annotation", label = FALSE, pt.size = 0.5) # + NoLegend()
# Plot_UMAP_Annotation

Plot_UMAP_Annotation <- DimPlot(seuratObject_Sample, reduction = "umap", group.by = "Annotation", label = FALSE, pt.size = 3) +
  scale_color_manual(values = color_CellType) +
  NoLegend() +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
        text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        aspect.ratio = 1)

Plot_UMAP_Annotation





# Plot_UMAP_ACT <- DimPlot(seuratObject_Sample, reduction = "umap", group.by = "Actual_Cell_Type", label = FALSE, pt.size = 0.5) # + NoLegend()
# Plot_UMAP_ACT

Plot_UMAP_ACT <- DimPlot(seuratObject_Sample, reduction = "umap", group.by = "Actual_Cell_Type", label = FALSE, pt.size = 3) +
  scale_color_manual(values = color_CellType) +
  NoLegend() +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
        text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        aspect.ratio = 1)

Plot_UMAP_ACT





## seuratObject_Ref
seuratObject_Ref@meta.data$`Actual_Cell_Type` <- seuratObject_Ref@meta.data$`Actual Cell Type`

# Plot_UMAP_ACT_Ref <- DimPlot(seuratObject_Ref, reduction = "umap", group.by = "Actual_Cell_Type", label = FALSE, pt.size = 0.5)
# Plot_UMAP_ACT_Ref



Plot_UMAP_ACT_Ref <- DimPlot(seuratObject_Ref, reduction = "umap", group.by = "Actual_Cell_Type", label = FALSE, pt.size = 3) +
  scale_color_manual(values = color_CellType) +
  NoLegend() +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2),
        text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 18),
        aspect.ratio = 1)

Plot_UMAP_ACT_Ref


Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)

# pdf(paste0(export_folder, "/", export_name, "_Graphical_Abstract.pdf"),
pdf(paste0(Name_time_wo_micro, "_Graphical_Abstract.pdf"),
    width = 16, height = 8)

print(Plot_UMAP_Unknown + Plot_UMAP_VICTORS)
# print(Plot_UMAP_VICTORS)
print(Plot_UMAP_Annotation)
print(Plot_UMAP_ACT)

print(Plot_UMAP_ACT_Ref)


dev.off()
