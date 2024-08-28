library(ggplot2)
library(cowplot)
library(dplyr)
library(rlang)

# 提取所有子图的 x 轴标签（细胞类型）
cell_types <- unique(unlist(lapply(plots1, function(plot) {
  x_var <- as_label(plot$mapping$x)  # 使用 as_label 提取变量名
  plot$data[[x_var]]
})))

# 创建编号与细胞类型的对照表
cell_type_mapping <- data.frame(
  Number = seq_along(cell_types),
  Cell_Type = cell_types
)

# 定義一個函數來處理每個 plot 的更新
update_plot <- function(plot, index) {
  plot <- plot + scale_x_discrete(labels = setNames(cell_type_mapping$Number, cell_types)) +
    theme(axis.text.x = element_text(angle = 0, size = 14),
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +  # 完全移除边距，减少空隙
    labs(title = gsub(" on Actual_Cell_Type", "", plot$labels$title)) +
    theme(plot.title = element_text(size = 14, face = "bold"))  # 放大标题字体

  if (index != 1) {
    plot <- plot + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }

  plot <- plot + theme(axis.title.x = element_blank())

  return(plot)
}

# 更新 plots1 和 plots2
plots1_1 <- c(plots[1:(length(plots)/2)])
plots1_1 <- lapply(seq_along(plots1_1), function(i) update_plot(plots1_1[[i]], i))
plots1 <- lapply(seq_along(plots1), function(i) update_plot(plots1[[i]], i))

plots2_1 <- c(plots[1:(length(plots)/2)])
plots2_1 <- lapply(seq_along(plots2_1), function(i) update_plot(plots2_1[[i]], i))
plots2 <- lapply(seq_along(plots2), function(i) update_plot(plots2[[i]], i))



# 合并并打印 combined_Percent_plot 和 combined_Percent_plot2
combined_plots1_1 <- plot_grid(
  plotlist = plots1_1,
  ncol = 7,
  align = 'hv',
  axis = "tb",
  rel_widths = rep(1, length(plots1_1)),
  rel_heights = rep(1, length(plots1_1))
)
print(combined_plots1_1)

combined_plots2_1 <- plot_grid(
  plotlist = plots2_1,
  ncol = 7,
  align = 'hv',
  axis = "tb",
  rel_widths = rep(1, length(plots2_1)),
  rel_heights = rep(1, length(plots2_1))
)
print(combined_plots2_1)



combined_Percent_plot <- plot_grid(
  plotlist = plots1,
  ncol = 7,
  align = 'hv',
  axis = "tb",
  rel_widths = rep(1, length(plots1)),
  rel_heights = rep(1, length(plots1))
)
print(combined_Percent_plot)

combined_Percent_plot2 <- plot_grid(
  plotlist = plots2,
  ncol = 7,
  align = 'hv',
  axis = "tb",
  rel_widths = rep(1, length(plots2)),
  rel_heights = rep(1, length(plots2))
)
print(combined_Percent_plot2)

# 输出编号与细胞类型的对照表
print(cell_type_mapping)
