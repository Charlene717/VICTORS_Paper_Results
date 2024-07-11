#### Check Markers ####

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


##### Load Package #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)

if(!require(cowplot)) install.packages("cowplot"); library(cowplot)

source("FUN_Plot_Beautify_UMAP_Box.R")
source("Set_plot_color.R")


#### Load Data ####
# Dataset <- "Pancreas"
# Load_Path <- "D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Pancreas/IntCT_scRNAseq_Ref"

# Dataset <- "GSE132044"
# Load_Path <- "D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/GSE132044/GSE132044_Read_All_Processed_Ref_OriQC"

Dataset <- "PMID37027478"
Load_Path <- "D:/Dropbox/##_GitHub/###_VUMC/CreateDataset/Input_Dataset/Lung/PMID37027478/PMID37027478_6paper_small_5000_PropSame_Flt"


file_list <- list.files(Load_Path, pattern = "\\.RData$", full.names = TRUE)


combined_plot_ncol <- 3
MarkerUMAP_Width <- 15
MarkerUMAP_height <- 10

# 初始化用於儲存UMAP圖的列表
umap_plots <- list()
bar_plots <- list()
percent_plots <- list()

# 循環讀取每個RData文件
for(file in file_list) {
  # 載入RData文件
  load(file)

  # seuratObject@meta.data$Actual_Cell_Type %>% unique()
  try({
    cells_to_keep <- rownames(seuratObject@meta.data)[!is.na(seuratObject@meta.data$Actual_Cell_Type)]
    seuratObject <- seuratObject[, cells_to_keep]

    seuratObject <- subset(seuratObject, subset = Actual_Cell_Type != "Unknown")
  })


  df <- data.frame(Actual_Cell_Type = as.character(seuratObject$Actual_Cell_Type),
                   seurat_clusters = as.character(seuratObject$seurat_clusters))

  if(Dataset == "Pancreas"){
    Name_Title <- seuratObject@misc[["BasicInfo"]][["DataID"]]
  }else if(Dataset == "GSE132044"){
    Name_Title <- seuratObject@misc[["BasicInfo"]][["Platform"]]
  }else{
    Name_Title <- seuratObject@misc[["BasicInfo"]][["DataID"]]
  }

  plots_Act <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "Actual_Cell_Type", Set_cluster_Title = Name_Title,
                                 Set_cluster2 = "Actual_Cell_Type",Set_cluster_Title2 = Name_Title,
                                 palette = "Set3", legend = FALSE, color_vector = color_CellType)

  print(plots_Act$UMAP_cluster + plots_Act$UMAP_label +
        plots_Act$Grouped_Barchart + plots_Act$Percent_Stacked_Barchart)

  umap_plots[[file]] <- plots_Act$UMAP_label
  bar_plots[[file]] <- plots_Act$Grouped_Barchart
  percent_plots[[file]] <- plots_Act$Percent_Stacked_Barchart

  # # 假設每個RData包含一個名為seuratObject的Seurat物件
  # # 生成UMAP圖，並用實際細胞類型標記
  # umap_plot <- DimPlot(seuratObject, reduction = "umap", group.by = "Actual_Cell_Type")
  #
  # # 添加UMAP圖到列表
  # umap_plots[[file]] <- umap_plot
}

## 合併所有UMAP圖成一個組圖
# library(gridExtra)
# combined_umap_plot <- do.call(grid.arrange, umap_plots)
# print(combined_umap_plot)

library(cowplot)
combined_plot_umap <- plot_grid(plotlist = umap_plots, ncol = 3)
print(combined_plot_umap)

combined_plot_bar <- plot_grid(plotlist = bar_plots, ncol = 3)
print(combined_plot_bar)

combined_plot_percent <- plot_grid(plotlist = percent_plots, ncol = 3)
print(combined_plot_percent)


#### Export ####
Name_ExportFolder <- "Export_CTUMAP_CellCount"

Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_Export <- paste0(Name_time_wo_micro, "_", Dataset)


if(Dataset == "Pancreas"){
  MarkerUMAP_Width <- 15
  MarkerUMAP_height <- 10
}else if(Dataset == "GSE132044"){
  MarkerUMAP_Width <- 10
  MarkerUMAP_height <- 15
}else if(Dataset == "PMID37027478"){
  MarkerUMAP_Width <- 10
  MarkerUMAP_height <- 8
}



if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_MarkerUMAP.pdf"),
    width = MarkerUMAP_Width, height = MarkerUMAP_height)
print(combined_plot_umap)
print(combined_plot_bar)
print(combined_plot_percent)
dev.off()

