### ChatGPT Record: https://chat.openai.com/share/105e6dcb-3092-4f95-9cad-e3bac193ae60

if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("reshape2")) install.packages("reshape2"); library(reshape2)
# if(!require("Seurat")) install.packages("Seurat"); library(Seurat)

try({source("Set_pbmc3k_plot.R")})


plot_box_PvStat <- function(data, cols, x = "Cell Type", x_label = "Cell Type", y_label = "Value",
                            legend_label = "Score Type", title = "Plot Title",
                            color_vector_Ref = NULL, score_suffix = " score") {

  # Convert data frame format to fit ggplot2
  df_melted <- melt(data[, c(x, cols)], id.vars = x)
  df_melted$variable <- gsub(paste0(score_suffix, "$|^pvalue of |^padj of "), "", df_melted$variable)
  # df_melted$variable <- gsub(" score$|^pvalue of |^padj of ", "", df_melted$variable)
  CellType.set <- df_melted$variable %>% unique()
  if (!is.null(color_vector_Ref) &&  all(CellType.set %in% names(color_vector_Ref))) {
    fill_colors <- scale_fill_manual(values = color_vector_Ref)
  } else {
    fill_colors <- NULL
  }

  p <- ggplot(df_melted, aes(x = !!sym(x), y = value, fill = variable)) +
    geom_boxplot() +
    labs(title = title, x = x_label, y = y_label, fill = legend_label) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14, angle = 0),
      title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      aspect.ratio = 1/3,
      panel.border = element_rect(color = "black", fill = NA, size = 1.5) # Ensure black border
    )

  if (!is.null(fill_colors)) {
    p <- p + fill_colors
  }

  return(p)
}


# ## Test Function
# # Extract meta_data
# meta_data <- seuratObject_Sample@meta.data
# meta_data_Ref <- seuratObject_Ref@meta.data
#
# # Grab the required column names
# score_cols <- grep("score$", colnames(meta_data), value=TRUE, ignore.case=TRUE)
# pvalue_cols <- grep("^pvalue.*score$", colnames(meta_data), value=TRUE, ignore.case=TRUE)
# padj_cols <- grep("^padj.*score$", colnames(meta_data), value=TRUE, ignore.case=TRUE)
#
# # Exclude pvalue and padj items from score_cols
# score_cols <- setdiff(score_cols, c(pvalue_cols, padj_cols))
#
#
# # Draw a box and whisker plot
# plot_box_CT1 <- plot_box_PvStat(meta_data, score_cols,x = "Cell Type", x_label = "Cell Type", y_label = "Module score", title = "Module Score", color_vector_Ref = color_vector_Ref)
# plot_box_CT1_Ref <- plot_box_PvStat(meta_data_Ref, score_cols,x = "Cell Type", x_label = "Cell Type", y_label = "Module score", title = "Reference Module Score", color_vector_Ref = color_vector_Ref)
# plot_box_CT2 <- plot_box_PvStat(meta_data, pvalue_cols,x = "Cell Type", x_label = "Cell Type", y_label = "pvalue", title = "P-value of Module Score", color_vector_Ref = color_vector_Ref)
# plot_box_CT3 <- plot_box_PvStat(meta_data, padj_cols,x = "Cell Type", x_label = "Cell Type", y_label = "padj", title = "Padj of Module Score", color_vector_Ref = color_vector_Ref)
#
# plot_box_CT1 / plot_box_CT1_Ref
# plot_box_CT2 / plot_box_CT3
#

# # Combining three box-and-whisker plots using grid.arrange
# # grid.arrange(p1, p2, p3, ncol=1)
# # gridExtra::grid.arrange(plot_box_CT1, plot_box_CT2, plot_box_CT3, ncol=3)


plot_box_PvStat_2 <- function(data, score_columns, cell_type_col = "Actual Cell Type", source_col = "Source", Title = "", color_vector_Ref = NULL) {
  require(ggplot2)
  require(dplyr)
  require(tidyr)

  # 將數據轉為長格式
  df_long <- data %>%
    select(all_of(c(cell_type_col, source_col, score_columns))) %>%
    gather(key = "CellType_Score", value = "Score", -c(cell_type_col, source_col))

  # 從CellType_Score列提取variable資訊
  df_long$variable <- gsub(" score$", "", df_long$CellType_Score)

  # 根據color_vector_Ref判定顏色
  CellType.set <- unique(df_long$variable)
  if (!is.null(color_vector_Ref) && all(CellType.set %in% names(color_vector_Ref))) {
    fill_colors <- scale_fill_manual(values = color_vector_Ref)
    fill_colors2 <- scale_color_manual(values = color_vector_Ref)
  } else {
    fill_colors <- NULL
    fill_colors2 <- NULL
  }

  # 使用ggplot繪圖
  p <- ggplot(df_long, aes(x = !!sym(cell_type_col), y = Score, color = variable, fill = variable)) +
    geom_boxplot() +
    facet_wrap(as.formula(paste("~", source_col)), scales = "free_x") +
    labs(title = Title, y = "Score", color = "Score Type", fill = "Score Type") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14, angle = 0),
      title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      aspect.ratio = 1/2,
      panel.border = element_rect(color = "black", fill = NA, size = 1.5) # Ensure black border
    )

  # 加入顏色調整
  if (!is.null(color_vector_Ref) && all(CellType.set %in% names(color_vector_Ref))) {
    p <- p + fill_colors + fill_colors2
  }

  print(p)
}

# ## Test Function
# # 提取結尾為 "score" 但開頭不是 "pvalue" 或 "padj" 的列
# metadata_all_C <- combined_obj@meta.data
# score_columns <- grep("^(?!pvalue|padj).*score$", names(metadata_all_C), value = TRUE, perl = TRUE)
#
# plot_box_PvStat_2(data = metadata_all_C, score_columns = score_columns,
#                   Title = "Module Score", color_vector_Ref = color_vector_Ref)
# plot_box_PvStat_2(data = metadata_all_C, score_columns = score_columns,cell_type_col = "Cell Type",
#                   Title = "Module Score", color_vector_Ref = color_vector_Ref)
