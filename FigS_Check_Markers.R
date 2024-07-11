## Ref: http://117.50.127.228/CellMarker/CellMarker_help.html

#### Check Markers ####

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


##### Load Package #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require(cowplot)) install.packages("cowplot"); library(cowplot)


#### Load Data ####
# Dataset <- "pbmc3k"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/pbmc3k/pbmc3k_Ori/Seurat_pbmc3k.RData")
#
# Dataset <- "PMID37027478_GSE122960"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Lung/PMID37027478/PMID37027478_6paper_small_5000/GSE122960_Sample_CellNum5000_Seed123_Processed.RData")


## GSE132044
# Dataset <- "GSE132044_10xV2"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Ori/GSE132044_10xV2.RData")
#
# Dataset <- "GSE132044_10xV2A"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Ori/GSE132044_10xV2A.RData")
#
# Dataset <- "GSE132044_10xV2B"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Ori/GSE132044_10xV2B.RData")
#
# Dataset <- "GSE132044_10xV3"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Ori/GSE132044_10xV3.RData")
#
# Dataset <- "GSE132044_CELSeq2"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Ori/GSE132044_CELSeq2.RData")
#
# Dataset <- "GSE132044_DropSeq"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Ori/GSE132044_DropSeq.RData")
#
# Dataset <- "GSE132044_inDrops"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Ori/GSE132044_inDrops.RData")
#
# Dataset <- "GSE132044_SeqWell"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Ori/GSE132044_SeqWell.RData")
#
# Dataset <- "GSE132044_Smartseq2"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Ori/GSE132044_Smartseq2.RData")

## Pancreas
# Dataset <- "Pancreas_Baron"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Pancreas/scRNAseq_Ref/BaronPancreasData.RData")
# Dataset <- "Pancreas_Muraro"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Pancreas/scRNAseq_Ref/MuraroPanData.RData")
# Dataset <- "Pancreas_Segerstolpe"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Pancreas/scRNAseq_Ref/SegerstolpePanData.RData")

# Dataset <- "Pancreas_Baron"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Pancreas/IntCT_scRNAseq_Ref/BaronPancreasData.RData")
# Dataset <- "Pancreas_Muraro"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Pancreas/IntCT_scRNAseq_Ref/MuraroPancreasData.RData")
# Dataset <- "Pancreas_Segerstolpe"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Pancreas/IntCT_scRNAseq_Ref/SegerstolpePancreasData.RData")

# Dataset <- "PMID37027478_GSE128033"
# load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Lung/PMID37027478/PMID37027478_GSE128033_small_5000_PropSame/GSE128033_Sample_CellNum5000_Seed123_Processed.RData")

Dataset <- "PMID37027478_GSE135893"
load("D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Lung/PMID37027478/PMID37027478_GSE135893_small_5000_PropSame/GSE135893_Sample_CellNum5000_Seed123_Processed.RData")

seuratObject@meta.data$Actual_Cell_Type %>% unique()

try({
  cells_to_keep <- rownames(seuratObject@meta.data)[!is.na(seuratObject@meta.data$Actual_Cell_Type)]
  seuratObject <- seuratObject[, cells_to_keep]
})

combined_plot_ncol <- 4
MarkerUMAP_Width <- 12
MarkerUMAP_height <- 10

## Define cell types and corresponding marker genes
if(grepl("GSE132044", Dataset)) {
  cell_types_markers <- list(
    "B cell" = c("MS4A1"), # "B cell" = c("CD19", "MS4A1"),
    "CD4+ T cell" = c("IL7R"), # "CD4+ T cell" = c("CD4", "IL7R"),
    "Cytotoxic T cell" = c("CD8A"), # "Cytotoxic T cell" = c("CD8A", "GZMB"),
    "Natural killer cell" = c("GNLY"), # "Natural killer cell" = c("GNLY", "NKG7"),
    "CD14+ monocyte" = c("CD14"), # "CD14+ monocyte" = c("CD14", "LYZ"),
    "CD16+ monocyte" = c("FCGR3A"), #    "CD16+ monocyte" = c("FCGR3A", "CD16"),
    "Dendritic cell" = c("CD1C"), # "Dendritic cell" = c("CD1C", "CD141"),
    "Megakaryocyte" = c("PPBP"), # "Megakaryocyte" = c("PPBP", "ITGA2B"),
    "Plasmacytoid dendritic cell" = c("TCF4") # "Plasmacytoid dendritic cell" = c("TCF4", "IL3RA")
  )
  # combined_plot_ncol <- 5
  # MarkerUMAP_Width <- 18
  # MarkerUMAP_height <- 6

  combined_plot_ncol <- 3
  MarkerUMAP_Width <- 12
  MarkerUMAP_height <- 10

}else if(grepl("PMID37027478", Dataset)){
  cell_types_markers <- list(
    "B cell" = c("MS4A1"), # "B cell" = c("CD19", "MS4A1"),
    "Endothelial cell" = c("PECAM1"),    # "Endothelial cell" = c("PECAM1", "VWF"),
    "Epithelial cell" = c("KRT18"), # "Epithelial cell" = c("EPCAM", "KRT18"),
    "Fibroblast" = c("COL1A1"), #     "Fibroblast" = c("COL1A1", "COL1A2"),
    "Lymphoid endothelial cell" = c("PDPN"), # "Lymphoid endothelial cell" = c("PECAM1", "LYVE1"),
    "Mast cell" = c("TPSAB1"), # "Mast cell" = c("TPSAB1", "KIT"),
    "Myeloid cell" = c("LYZ"), # "Myeloid cell" = c("LYZ", "CD68"),
    "Natural killer T cell" = c("TRAC"), # "Natural killer T cell" = c("XCL1", "XCL2"), # "CD3E"
    "Proliferative cell" = c("MKI67"), #  "Proliferative cell" = c("MKI67", "PCNA"),
    "Smooth muscle Cell" = c("ACTA2"), # "Smooth muscle Cell" = c("ACTA2", "MYH11"),
    "Type I pneumocyte" = c("EMP2"), # "Type I pneumocyte" = c("AGER", "PDPN"),
    "Type II pneumocyte" = c("SFTPC") # "Type II pneumocyte" = c("SFTPC", "ABCA3")
  )

  combined_plot_ncol <- 3
  MarkerUMAP_Width <- 10
  MarkerUMAP_height <- 12

}else if(grepl("pbmc3k", Dataset)){
  cell_types_markers <- list(
    "Memory CD4 T" = c("IL7R", "CCR7"),
    "B" = c("CD19", "MS4A1"),
    "CD14+ Mono" = c("CD14", "LYZ"),
    "NK" = c("GNLY", "NKG7"),
    "CD8 T" = c("CD8A", "CD8B"),
    "Naive CD4 T" = c("CCR7", "SELL"),
    "FCGR3A+ Mono" = c("FCGR3A", "CD16"),
    "DC" = c("CD1C", "CD141"),
    "Platelet" = c("PPBP", "PF4")
  )


}else if(grepl("Pancreas", Dataset)){
  cell_types_markers <- list(
    "Acinar cell" = c("PRSS1"),  # "Acinar cell" = c("PRSS1", "CPA1"),
    "Alpha cell" = c("GCG"), # "Alpha cell" = c("GCG", "ARX"),
    "Beta cell" = c("INS"), # "Beta cell" = c("INS", "PDX1"),
    "Delta cell" = c("SST"), # "Delta cell" = c("SST", "HHEX"),
    "Ductal cell" = c("KRT19"), # "Ductal cell" = c("KRT19", "SOX9"),
    "Endothelial cell" = c("PECAM1"),  # "Endothelial cell" = c("PECAM1", "VWF"),
    # "Stellate cell" = c("COL1A1", "ACTA2"),
    # "Gamma cell" = c("PPY"),
    # "Macrophage cell" = c("CD68", "CD14"),
    # "Schwann cell" = c("S100B", "MBP"),
    # "Mast cell" = c("TPSAB1", "KIT"),
    # "T cell" = c("CD3D", "CD3E"),
    "Epsilon cell" = c("GHRL")
  )
  combined_plot_ncol <- 7
  MarkerUMAP_Width <- 22
  MarkerUMAP_height <- 3

}



## Create a list to store UMAP
umap_plots <- list()

## Plot UMAP for each cell type and corresponding marker gene
for(cell_type in names(cell_types_markers)) {
  for(marker in cell_types_markers[[cell_type]]) {
    # Use tryCatch to catch and handle possible errors
    p <- tryCatch({
      FeaturePlot(seuratObject, features = marker, cols = c("lightgrey", "blue")) +
        ggtitle(paste(marker, "(", cell_type, ")", sep = ""))
    }, error = function(e) {
      # If an error occurs, return NULL and print the error message
      message("Error plotting ", marker, " in ", cell_type, ": ", e$message)
      NULL
    })

    # If p is not NULL, add to umap_plots list
    if(!is.null(p)) {
      umap_plots[[paste(cell_type, marker, sep = "_")]] <- p
    }
  }
}


## Combine all UMAP
combined_plot <- plot_grid(plotlist = umap_plots, ncol = combined_plot_ncol)
print(combined_plot)

##
p1 <- DimPlot(seuratObject, label = TRUE, repel = TRUE, group.by = "seurat_clusters") + NoLegend()
p2 <- DimPlot(seuratObject, label = TRUE, repel = TRUE, group.by = "Actual_Cell_Type") + NoLegend()
p1 + p2
combined_plot2 <- plot_grid(plotlist = list(combined_plot, p1 + p2), ncol = 2)

## Export
Name_ExportFolder <- "Export_Check_Markers"

Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_Export <- paste0(Name_time_wo_micro, "_", Dataset)

if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MarkerUMAP.pdf"),
    width = MarkerUMAP_Width, height = MarkerUMAP_height)
print(combined_plot)
dev.off()
