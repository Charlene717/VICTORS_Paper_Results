if(!require("tidyr")) install.packages("tidyr"); library(tidyr)
if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if(!require("dplyr")) install.packages("dplyr"); library(dplyr)
if(!require("patchwork")) install.packages("patchwork"); library(patchwork)

## Set Para
Set_boxplot_fill <- FALSE # Set_boxplot_fill <- TRUE


## Functions of plot
my_ggplot_theme <- function() {
  theme_bw() +
    theme(
      text = element_text(size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
      axis.title = element_text(size = 24),
      axis.text.y = element_text(size = 22),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 20),
      panel.border = element_rect(colour = "black", fill = NA, size = 1.2)  # 加粗外围边框
    )
}

if(Set_boxplot_fill){
  create_metric_plot <- function(data, metric, figure_note, color_method, x_col = "Actual_Cell_Type", set_ylimits = c(0, 1)) {
    ggplot(data %>% filter(Metric == metric),
           aes(x = .data[[x_col]], y = Value, fill = Method)) +
      geom_boxplot() +
      # scale_y_continuous(limits = c(0, 1)) +
      my_ggplot_theme() +
      labs(title = figure_note, # paste0(figure_note, " ", metric, " Across Actual Cell Types"),
           x = x_col, y = metric, fill = "Method") +
      scale_fill_manual(values = color_method)
  }
}else{
  create_metric_plot <- function(data, metric, figure_note, color_method, x_col = "Actual_Cell_Type", set_ylimits = c(0, 1)) {
    # 确保数据类型正确
    data$Value <- as.numeric(data$Value)
    data[[x_col]] <- factor(data[[x_col]])
    # data$Method <- factor(data$Method, levels = names(color_method))

    # 创建盒形图
    plt <- ggplot(data %>% filter(Metric == metric),
                  aes(x = .data[[x_col]], y = Value)) +
      geom_boxplot(aes(color = Method, fill = Method),alpha = 0.3,
                   outlier.shape = 21, outlier.colour = NA, size = 0.8) +  # 明确边框颜色
      # scale_y_continuous(limits = c(0, 1)) +
      scale_fill_manual(values = color_method) +
      scale_color_manual(values = color_method) +
      theme_classic() +
      theme(
        panel.background = element_rect(fill = NA, colour = "black", size = 1.3),
        panel.grid.major.x = element_blank()
      ) +
      my_ggplot_theme() +
      labs(title = figure_note,
           x = x_col, y = metric, color = "Method")

    # Add outlier
    plt <- plt + geom_point(aes(color = Method), position = position_dodge(width = 0.75), size = 1.2)
    # plt <- plt + stat_summary(fun = mean, geom = "text", aes(label = round(..y.., digits = 2)),
    #                           vjust = -1.5, size = 6, color = "black")

    return(plt)
  }

}


## Set Color
color_Method <- list(
  "VICTORS" =  "#32383b", # "#044a63", # "#2c7856",

  "singleR" = "#2862bf",   # "#884db8" , # "#b853b1", # "#2c7856",
  "singleR_VICTOR" = "#32383b", # "#2c7856",
  # "singleR_VICTORS" = "#c4126a", # "#2c7856",

  "scPred" = "#368a5b", # "#13c274",
  "scPred_VICTOR" = "#32383b", # "#2c7856",

  "SCINA" = "#e0385d", # "#8f174d",
  "SCINA_VICTOR" = "#32383b", # "#2c7856",

  "scmap" = "#ad772f", # "#3e75bd" # "#696866" # "#d94185"
  "scmap_VICTOR" = "#32383b", # "#2c7856",

  "CHETAH" = "#7ccf7f", # "#3e75bd" # "#696866" # "#d94185"
  "CHETAH_VICTOR" = "#32383b", # "#2c7856",

  "scClassify" = "#8d389c", # "#3e75bd" # "#696866" # "#d94185"
  "scClassify_VICTOR" = "#32383b", # "#2c7856",

  "Seurat" = "#42dbd9", # "#3e75bd" # "#696866" # "#d94185"
  "Seurat_VICTOR" = "#32383b" # "#2c7856"
)

# source("PlotFun_SetColor.R")
# if (!exists("color_Method")) {color_Method <- setNames(character(0), character(0))}
# color_Method <- update_color(c(data_10x$Method,"VICTORS"), color_Method)

## create plots list
create_plots <- function(data, methods, metric, color_method, set_x_col = "Actual_Cell_Type",
                         set_ylimits = NULL) { # set_ylimits = NULL # set_ylimits = c(0, 1)
  # 初始化一个空的列表来存储生成的图
  plots <- list()

  for (method in methods) {
    # 为每个方法生成过滤后的数据子集
    filtered_data <- data %>%
      filter(grepl(paste0("^", method), Method)) %>%
      mutate(Method = gsub(".*_VICTORS$", "VICTORS", Method))

    # 生成图表
    plot <- create_metric_plot(filtered_data, metric, paste0("Annotation by ", method), color_method, x_col = set_x_col, set_ylimits = c(0, 1))
    if(!is.null(set_ylimits)){
      plot <- plot + scale_y_continuous(limits = set_ylimits) # plot <- plot + scale_y_continuous(limits = c(0, 1))
    }

    plots[[method]] <- plot
  }

  return(plots)
}

# ## Test function
# methods <- c("singleR", "scmap", "SCINA", "scPred")
# plots <- create_plots(data_10x, methods, "Accuracy", color_Method)
#
# combined_plots <- wrap_plots(plots, ncol = 2)
# combined_plots

################################
#### Create Combind Plot ####
create_combind_plot <- function(plots, ncol, title_data_10x, legend_data, color_legend, legend_title = "Method") {
  # 检查并更新颜色定义
  if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
  if(!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
  if(!require("patchwork")) install.packages("patchwork"); library(patchwork)
  if(!require("cowplot")) install.packages("cowplot"); library(cowplot)

  source("PlotFun_SetColor.R")
  if (!exists("color_legend")) {color_legend <- setNames(character(0), character(0))}
  color_legend <- update_color(legend_data, color_legend)

  # 计算行数
  nrow <- ceiling(length(plots) / ncol)

  # 调整子图的轴标题和轴文本的可见性
  for (i in seq_along(plots)) {
    current_row <- ceiling(i / ncol)
    current_col <- ifelse(i %% ncol == 0, ncol, i %% ncol)

    # 应用主题元素
    plots[[i]] <- plots[[i]] + theme(text = element_text(size = 20), legend.position = "none",
                                     axis.title.x = element_text(face = "bold", size = 24),
                                     axis.title.y = element_text(face = "bold", size = 24))
    if (current_col != 1) {
      # plots[[i]] <- plots[[i]] + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
      #                                  axis.ticks.y = element_blank())
      plots[[i]] <- plots[[i]] + theme(axis.title.y = element_blank())
    }
    if (current_row != nrow) {
      plots[[i]] <- plots[[i]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                                       axis.ticks.x = element_blank())
    }
  }


  # Create a dummy data frame for the legend with a single point
  legend_data <- data.frame(Method = names(color_legend), x = 1, y = 1)
  custom_order <- c("singleR",  "scmap", "SCINA", "scPred",
                    "CHETAH","scClassify","Seurat","VICTOR" )
  legend_data$Method <- factor(legend_data$Method, levels = custom_order)

  # Create a ggplot object for the legend without actual points
  legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Method)) +
    geom_point(shape = 15, size = 0, show.legend = TRUE) +  # Use shape 15 for squares
    scale_color_manual(values = color_legend) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 22)) +
    guides(color = guide_legend(title = legend_title, override.aes = list(shape = 15,size =6)))

  # 创建标题标签
  title_label <- ggdraw() + draw_plot_label(title_data_10x, size = 15, x = 0.5, hjust = 0.5)

  # 组合所有图元素
  final_plot <- plot_grid(title_label, wrap_plots(plots, ncol = ncol), legend_plot, ncol = 1, rel_heights = c(0.1, 1, 0.1))

  return(final_plot)
}

## Set Color
# Define the colors for the methods
color_legend <- c(
  "VICTOR" = "#32383b", # "#044a63"
  "singleR_VICTOR" = "#32383b", # "#044a63"
  "scPred_VICTOR" = "#32383b", # "#044a63"
  "SCINA_VICTOR" = "#32383b", # "#044a63"
  "scmap_VICTOR" = "#32383b", # "#044a63"
  "CHETAH_VICTOR" = "#32383b", # "#044a63"
  "scClassify_VICTOR" = "#32383b", # "#044a63"
  "Seurat_VICTOR" = "#32383b", # "#044a63"

  "singleR" = "#2862bf",   # "#884db8"
  "scPred" = "#368a5b",
  "SCINA" = "#e0385d",
  "scmap" = "#ad772f",
  "CHETAH" = "#7ccf7f",
  "scClassify" = "#8d389c",
  "Seurat" = "#42dbd9"
)

# source("PlotFun_SetColor.R")
# if (!exists("color_legend")) {color_legend <- setNames(character(0), character(0))}
# color_legend <- update_color(c(gsub(".*_VICTORS$", "VICTORS", as.character(data_10x$Method)),"VICTORS"), color_legend)


# ## Example usage:
# legend_set <- c(gsub(".*_VICTORS$", "VICTORS", as.character(data_10x$Method))) %>% unique()
# final_plot <- create_combind_plot(plots, 2, Title_data_10x, legend_set, color_legend)
# print(final_plot)

## create_and_combine_metric_plots
create_and_combine_metric_plots <- function(data, methods, figure_note, metric, platform_type, ncol, legend_set, color_legend, set_x_col = "Actual_Cell_Type", set_ylimits = c(0, 1)) {
  title <- paste0("Evaluation on ", figure_note, ": ", metric, " Across Actual Cell Types - ", platform_type)
  plots <- create_plots(data, methods, metric, color_legend, set_x_col = set_x_col, set_ylimits = set_ylimits)
  final_plot <- create_combind_plot(plots, ncol, title, legend_set, color_legend)
  return(final_plot)
}

# ## Test function
# # 初始化方法和图例设置
# methods <- c("singleR", "scmap", "SCINA", "scPred")
# legend_set <- unique(gsub(".*_VICTORS$", "VICTORS", as.character(long_data$Method)))
#
# # 示例：为不同的数据集和指标创建组合图
# final_plots_data_10x_Accuracy <- create_and_combine_metric_plots(data_10x, methods, Figure_Note, "Accuracy", "10x Platform", 2, legend_set, color_legend)
# print(final_plots_data_10x_Accuracy)


