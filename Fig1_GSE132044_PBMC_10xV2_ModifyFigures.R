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

# 更新 plots 的函数
update_plots <- function(plots) {
  lapply(seq_along(plots), function(i) update_plot(plots[[i]], i))
}

# 将 plots 分为两组并更新
half_length <- length(plots) / 2
plots1_1 <- update_plots(plots[1:half_length])
plots1 <- update_plots(plots1)

plots2_1 <- update_plots(plots[1:half_length])
plots2 <- update_plots(plots2)

# 合并并打印 combined plots
combine_and_print <- function(plots, ncol = 7) {
  combined_plot <- plot_grid(
    plotlist = plots,
    ncol = ncol,
    align = 'hv',
    axis = "tb",
    rel_widths = rep(1, length(plots)),
    rel_heights = rep(1, length(plots))
  )
  print(combined_plot)
  return(combined_plot)
}

combined_plots1_1 <- combine_and_print(plots1_1)
combined_plots2_1 <- combine_and_print(plots2_1)
combined_Percent_plot <- combine_and_print(plots1)
combined_Percent_plot2 <- combine_and_print(plots2)

# 打印组合的 plots
print(combined_plots1_1 / combined_Percent_plot)
print(combined_plots2_1 / combined_Percent_plot2)


# 输出编号与细胞类型的对照表
print(cell_type_mapping)
