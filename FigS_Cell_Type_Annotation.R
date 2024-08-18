#### Check Markers ####

##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)


##### Load Package #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require(cowplot)) install.packages("cowplot"); library(cowplot)

# source("PlotFun_Beautify_UMAP_Box.R")
source("FUN_Plot_Beautify_UMAP_Box.R")
source("Set_plot_color.R")


#### Load Dataset ####
Dataset <- "GSE132044_PBMC_MislabelB" # "GSE132044_PBMC_MislabelB" # GSE132044_PBMC_MislabelNone"

if(Dataset == "GSE132044_PBMC_MislabelB"){
  # load("D:/Dropbox/###_VUMC/##_Research/VICTORS/20231229_Figures/PBMC_GSE132044/Export_GSE132044_MislabelB cell/20231212125506KYHDNV_Multi/20231212125506KYHDNV.RData")
  load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/#_Export_20240717/Export_Fig1_2024071601VAI_GSE132044_MisLabelB_Ref10xV2A/Fig1_2024071601VAI_GSE132044_MisLabelB_Fig1_Accuracy.RData")

}else if(Dataset == "GSE132044_PBMC_MislabelNone"){
  # load("D:/Dropbox/###_VUMC/##_Research/VICTORS/20231229_Figures/PBMC_GSE132044/Export_GSE132044_MislabelNone/20231221065523SVLJPQ_Multi/20231221065523SVLJPQ.RData")
  load("D:/Dropbox/##_GitHub/###_VUMC/VICTORS_Paper_Results/#_Export_20240717/20240712095102QILHSN_MislabelNone_Qry_10xV2_Ref_10xV2A/20240712095102QILHSN.RData")

}

seuratObject <- seuratObject_Sample
seuratObject@meta.data$`Actual Cell Type` %>% unique()

#### Set parameter ####
Name_time_wo_micro <- substr(gsub("[- :]", "", as.character(Sys.time())), 1, 14)
Name_FileID <- paste0(Name_time_wo_micro, paste0(sample(LETTERS, 3), collapse = ""))

Name_Export <- paste0(Name_FileID, "_", Dataset)

Name_ExportFolder <- paste0("Export_Annotation_",Name_Export)

if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder


#### Extract Annotation ####
seuratObject@meta.data$singleR <- seuratObject@meta.data$label_singleR
seuratObject@meta.data$scmap <- seuratObject@meta.data$label_scmap
seuratObject@meta.data$SCINA <- seuratObject@meta.data$label_SCINA
seuratObject@meta.data$scPred <- seuratObject@meta.data$label_scPred
seuratObject@meta.data$CHETAH <- seuratObject@meta.data$label_CHETAH
seuratObject@meta.data$scClassify <- seuratObject@meta.data$label_scClassify
seuratObject@meta.data$Seurat <- seuratObject@meta.data$label_Seurat


colnames(seuratObject@meta.data) <- gsub(" ", "_",colnames(seuratObject@meta.data))

df <- data.frame(`Actual_Cell_Type` = as.character(seuratObject$`Actual_Cell_Type`),
                 singleR = as.character(seuratObject$`singleR`),
                 scmap = as.character(seuratObject$`scmap`),
                 SCINA = as.character(seuratObject$`SCINA`),
                 scPred = as.character(seuratObject$`scPred`),
                 CHETAH = as.character(seuratObject$`CHETAH`),
                 scClassify = as.character(seuratObject$`scClassify`),
                 Seurat = as.character(seuratObject$`Seurat`),
                 `seurat_clusters` = as.character(seuratObject$`seurat_clusters`))


plots_Anno_CT_count1 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "singleR", Set_cluster_Title = "singleR",
                                         Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_CellType_Ref)
print(plots_Anno_CT_count1$UMAP_label2 + plots_Anno_CT_count1$UMAP_label + plots_Anno_CT_count1$Grouped_Barchart + plots_Anno_CT_count1$Percent_Stacked_Barchart)

plots_Anno_CT_count2 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "scmap", Set_cluster_Title = "scmap",
                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_CellType_Ref)
print(plots_Anno_CT_count2$UMAP_label2 + plots_Anno_CT_count2$UMAP_label + plots_Anno_CT_count2$Grouped_Barchart + plots_Anno_CT_count2$Percent_Stacked_Barchart)

plots_Anno_CT_count3 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "SCINA", Set_cluster_Title = "SCINA",
                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_CellType_Ref)
print(plots_Anno_CT_count3$UMAP_label2 + plots_Anno_CT_count3$UMAP_label + plots_Anno_CT_count3$Grouped_Barchart + plots_Anno_CT_count3$Percent_Stacked_Barchart)

plots_Anno_CT_count4 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "scPred", Set_cluster_Title = "scPred",
                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_CellType_Ref)
print(plots_Anno_CT_count4$UMAP_label2 + plots_Anno_CT_count4$UMAP_label + plots_Anno_CT_count4$Grouped_Barchart + plots_Anno_CT_count4$Percent_Stacked_Barchart)

plots_Anno_CT_count5 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "CHETAH", Set_cluster_Title = "CHETAH",
                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_CellType_Ref)
print(plots_Anno_CT_count5$UMAP_label2 + plots_Anno_CT_count5$UMAP_label + plots_Anno_CT_count5$Grouped_Barchart + plots_Anno_CT_count5$Percent_Stacked_Barchart)

plots_Anno_CT_count6 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "scClassify", Set_cluster_Title = "scClassify",
                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_CellType_Ref)
print(plots_Anno_CT_count6$UMAP_label2 + plots_Anno_CT_count6$UMAP_label + plots_Anno_CT_count6$Grouped_Barchart + plots_Anno_CT_count6$Percent_Stacked_Barchart)

plots_Anno_CT_count7 <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "Seurat", Set_cluster_Title = "Seurat",
                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_CellType_Ref)
print(plots_Anno_CT_count7$UMAP_label2 + plots_Anno_CT_count7$UMAP_label + plots_Anno_CT_count7$Grouped_Barchart + plots_Anno_CT_count7$Percent_Stacked_Barchart)


## Combine all plots
plots_UMAP.lt <- list()
plots_UMAP.lt[["UMAP1"]] <- plots_Anno_CT_count1$UMAP_label
plots_UMAP.lt[["UMAP2"]] <- plots_Anno_CT_count2$UMAP_label
plots_UMAP.lt[["UMAP3"]] <- plots_Anno_CT_count3$UMAP_label
plots_UMAP.lt[["UMAP4"]] <- plots_Anno_CT_count4$UMAP_label
plots_UMAP.lt[["UMAP5"]] <- plots_Anno_CT_count5$UMAP_label
plots_UMAP.lt[["UMAP6"]] <- plots_Anno_CT_count6$UMAP_label
plots_UMAP.lt[["UMAP7"]] <- plots_Anno_CT_count7$UMAP_label

plots_Percent.lt <- list()
plots_Percent.lt[["Percent1"]] <- plots_Anno_CT_count1$Percent_Stacked_Barchart+ theme(legend.position = "none")
plots_Percent.lt[["Percent2"]] <- plots_Anno_CT_count2$Percent_Stacked_Barchart+ theme(legend.position = "none")
plots_Percent.lt[["Percent3"]] <- plots_Anno_CT_count3$Percent_Stacked_Barchart+ theme(legend.position = "none")
plots_Percent.lt[["Percent4"]] <- plots_Anno_CT_count4$Percent_Stacked_Barchart+ theme(legend.position = "none")
plots_Percent.lt[["Percent5"]] <- plots_Anno_CT_count5$Percent_Stacked_Barchart+ theme(legend.position = "none")
plots_Percent.lt[["Percent6"]] <- plots_Anno_CT_count6$Percent_Stacked_Barchart+ theme(legend.position = "none")
plots_Percent.lt[["Percent7"]] <- plots_Anno_CT_count7$Percent_Stacked_Barchart+ theme(legend.position = "none")



combined_UMAP_plot <- plot_grid(plotlist = plots_UMAP.lt, ncol = 7)
print(combined_UMAP_plot)

combined_Percent_plot <- plot_grid(plotlist = plots_Percent.lt, ncol = 7)
print(combined_Percent_plot)


## 精簡整合圖
# 对 UMAP 图进行主题调整
plots_UMAP.lt <- lapply(seq_along(plots_UMAP.lt), function(i) {
  plot <- plots_UMAP.lt[[i]]

  # 对于除第一张图外的图，移除 y 轴标题
  if (i != 1) {
    plot <- plot + theme(axis.title.y = element_blank())
  }

  # 对所有图移除 x 轴标题和刻度，移除 y 轴刻度
  plot <- plot + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()  # 确保移除 y 轴刻度
  )

  return(plot)
})

# 合并 UMAP 图，确保第一张图和其他图大小一致
combined_UMAP_plot <- plot_grid(plotlist = plots_UMAP.lt, ncol = 7, align = 'hv')
print(combined_UMAP_plot)


## 重整 combined_Percent_plot
library(ggplot2)
library(cowplot)
library(dplyr)
library(rlang)

# 提取所有子图的 x 轴标签（细胞类型）
cell_types <- unique(unlist(lapply(plots_Percent.lt, function(plot) {
  x_var <- as_label(plot$mapping$x)  # 使用 as_label 提取变量名
  plot$data[[x_var]]
})))

# 创建编号与细胞类型的对照表
cell_type_mapping <- data.frame(
  Number = seq_along(cell_types),
  Cell_Type = cell_types
)

# 映射细胞类型到编号，并更新每个子图的 x 轴标签
plots_Percent.lt <- lapply(seq_along(plots_Percent.lt), function(i) {
  plot <- plots_Percent.lt[[i]]

  # 修改 x 轴标签为编号
  plot <- plot + scale_x_discrete(labels = setNames(cell_type_mapping$Number, cell_types)) +
    theme(axis.text.x = element_text(angle = 0, size = 14),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))  # 完全移除边距，减少空隙

  # 修改标题，去掉 "on Actual_Cell_Type"，并放大标题字体
  plot <- plot + labs(title = gsub(" on Actual_Cell_Type", "", plot$labels$title)) +
    theme(plot.title = element_text(size = 14, face = "bold"))  # 放大标题字体

  # 对于除第一张图外的图，移除 y 轴标题、刻度和标签
  if (i != 1) {
    plot <- plot + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }

  # 仅移除所有图的 x 轴标题
  plot <- plot + theme(axis.title.x = element_blank())

  return(plot)
})

# 合并 Percent 图，调整空隙，去除子图之间的空隙
combined_Percent_plot <- plot_grid(
  plotlist = plots_Percent.lt,
  ncol = 7,
  align = 'hv',
  axis = "tb",
  rel_widths = rep(1, length(plots_Percent.lt)),
  rel_heights = rep(1, length(plots_Percent.lt))
)
print(combined_Percent_plot)

# 输出编号与细胞类型的对照表
print(cell_type_mapping)

#### Create Figure legend ####
color_CellType_Ref <- list(
  "B cell" = "#ff7f0e",             # orange
  "CD14+ monocyte" = "#2ca02c",    # green
  "CD16+ monocyte" = "#e377c2",  # pink
  "CD4+ T cell" = "#1f77b4",  # blue
  "Cytotoxic T cell" = "#9467bd",         # purple
  "Dendritic cell" = "#bcbd22",            # olive
  "Megakaryocyte" = "#8c564b",   # brown
  "Natural killer cell" = "#d62728",            # red
  "Plasmacytoid dendritic cell" = "#17becf",       # light blue
  # "pruned" = "#f0ad8d",  #"#4a1919",       # dark grey
  "Unassign" = "#84e8d4"
)

library(ggplot2)
library(cowplot)

# 创建一个数据框来存储 cell_type_mapping 和颜色信息
legend_data <- data.frame(
  Cell_Type = names(color_CellType_Ref),
  Number = seq_along(color_CellType_Ref),
  Color = unlist(color_CellType_Ref)
)

# 创建一个显示编号和细胞类型的标签
legend_data$Label <- paste0(legend_data$Number, ": ", legend_data$Cell_Type)

# 确保按照 Number 的顺序进行排列
legend_data <- legend_data[order(legend_data$Number), ]

# 创建一个虚拟的 ggplot 对象来生成图例
legend_plot <- ggplot(legend_data, aes(x = factor(1), y = Number, fill = Label)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values = legend_data$Color, labels = legend_data$Label) +
  theme_void() +  # 移除所有非必要元素
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),  # 移除图例标题
    legend.text = element_text(size = 14),  # 设置图例文本大小
    legend.background = element_rect(fill = "white", color = NA),  # 图例背景设置为白色
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0.5, 'cm')  # 设置图例项之间的间距
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))  # 将图例设置为两列显示，并按照顺序排列

# 使用 cowplot 显示只有图例的图形
legend_only_plot <- plot_grid(legend_plot, ncol = 1)
print(legend_only_plot)

#### Export ####
pdf(paste0(Name_ExportFolder, "/", Name_Export, "_Annotation.pdf"),
    width = 15, height = 7)
print(combined_UMAP_plot)
print(combined_Percent_plot)
print(legend_only_plot)
dev.off()


# ## Count Cell type on cluster
# plots_CT_Clt_count <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "Actual_Cell_Type", Set_cluster_Title = "Actual_Cell_Type",
#                                         palette = "Set3", legend = FALSE, color_vector = color_vector)
# print(plots_CT_Clt_count$UMAP_cluster + plots_CT_Clt_count$UMAP_label + plots_CT_Clt_count$Grouped_Barchart + plots_CT_Clt_count$Percent_Stacked_Barchart)
#
# plots_Anno_CT_count <- Fun_Plot_UMAP_Bar(df, seuratObject, Set_cluster = "Actual_Cell_Type", Set_cluster_Title = "Actual_Cell_Type",
#                                          Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_vector)
# print(plots_Anno_CT_count$UMAP_label2 + plots_Anno_CT_count$UMAP_label + plots_Anno_CT_count$Grouped_Barchart + plots_Anno_CT_count$Percent_Stacked_Barchart)
#
#
# ## Export
# pdf(paste0(Name_ExportFolder,"/",Name_Export,"_CellTypeCount.pdf"),
#     width = 12, height = 12)
#
# print(plots_CT_Clt_count$UMAP_cluster + plots_CT_Clt_count$UMAP_label + plots_CT_Clt_count$Grouped_Barchart + plots_CT_Clt_count$Percent_Stacked_Barchart)
# print(plots_Anno_CT_count$UMAP_label2 + plots_Anno_CT_count$UMAP_label + plots_Anno_CT_count$Grouped_Barchart + plots_Anno_CT_count$Percent_Stacked_Barchart)
# print(plots_Anno_CT_count4$UMAP_label2 + plots_Anno_CT_count4$UMAP_label + plots_Anno_CT_count4$Grouped_Barchart + plots_Anno_CT_count4$Percent_Stacked_Barchart)
#
# dev.off()
#
#
#
# ###################
# #### singleR_Pruned ####
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if(!require("SingleR")) BiocManager::install("SingleR"); library(SingleR)
#
# Num_Col <- seuratObject@meta.data$`Actual_Cell_Type` %>% unique() %>% length()
# plots_DeltaDist <- plotDeltaDistribution(SingleR.lt, size = 2,ncol = Num_Col)
# print(plots_DeltaDist)
#
# seuratObject@meta.data$singleR_Pruned <- SingleR.lt@listData[["pruned.labels"]]
# seuratObject@meta.data$singleR_Pruned[is.na(seuratObject@meta.data$singleR_Pruned)] <- "pruned"
# df_singleR_Pruned <- data.frame(`Actual_Cell_Type` = as.character(seuratObject$`Actual_Cell_Type`),
#                                 singleR_Pruned = as.character(seuratObject$`singleR_Pruned`),
#                                 `seurat_clusters` = as.character(seuratObject$`seurat_clusters`))
#
# plots_Anno_CT_count1 <- Fun_Plot_UMAP_Bar(df_singleR_Pruned, seuratObject, Set_cluster = "singleR_Pruned", Set_cluster_Title = "singleR_Pruned",
#                                           Set_cluster2 = "Actual_Cell_Type", Set_cluster_Title2= "Actual_Cell_Type", palette = "Set3", legend = FALSE, color_vector = color_vector)
# print(plots_Anno_CT_count1$UMAP_label2 + plots_Anno_CT_count1$UMAP_label + plots_Anno_CT_count1$Grouped_Barchart + plots_Anno_CT_count1$Percent_Stacked_Barchart)
#
#
# ## Export
# pdf(paste0(Name_ExportFolder, "/", Name_Export, "_singleR_Pruned1.pdf"),
#     width = 15, height = 4)
# print(plots_DeltaDist)
# dev.off()
#
# pdf(paste0(Name_ExportFolder, "/", Name_Export, "_singleR_Pruned2.pdf"),
#     width = 12, height = 12)
# print(plots_Anno_CT_count1$UMAP_label2 + plots_Anno_CT_count1$UMAP_label + plots_Anno_CT_count1$Grouped_Barchart + plots_Anno_CT_count1$Percent_Stacked_Barchart)
#
# dev.off()
