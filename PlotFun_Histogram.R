if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

# Plot histogram with counts or proportions annotated
plot_histogram <- function(data, column, column2, type = "count", aspect_ratio = 1, Note_Title = "", Legend_Title = NULL, color_vector = NULL, position_type = "dodge") {

  all_combinations <- expand.grid(column_value = unique(data[[column]]), column2_value = unique(data[[column2]]))

  aggregated_data_count <- data %>%
    group_by(!!sym(column), !!sym(column2)) %>%
    summarise(Count = n()) %>%
    ungroup()

  aggregated_data_prop <- data %>%
    group_by(!!sym(column)) %>%
    mutate(Total = n()) %>%
    group_by(!!sym(column), !!sym(column2), Total) %>%
    summarise(Proportion = n() / unique(Total)) %>%
    ungroup()

  if (type == "count") {
    merged_data <- left_join(all_combinations, aggregated_data_count, by = c(column_value = column, column2_value = column2))
    plot <- ggplot(merged_data, aes_string(x = "column_value", y = "Count", fill = "column2_value")) +
      geom_bar(stat = "identity", position = position_type, width = 0.7) +
      (if (position_type == "dodge") geom_text(aes_string(label = "Count"), position = position_dodge(0.7), vjust = -0.5, angle = 45) else geom_text(aes_string(label = "Count"), position = position_stack(vjust = 0.5), angle = 45)) +
      labs(y = "Count", x = column, fill = column2)


  } else {
    merged_data <- left_join(all_combinations, aggregated_data_prop, by = c(column_value = column, column2_value = column2))
    plot <- ggplot(merged_data, aes_string(x = "column_value", y = "Proportion", fill = "column2_value")) +
      geom_bar(stat = "identity", position = "fill", width = 0.7) +
      geom_text(aes_string(label = "sprintf('%.1f%%', 100*Proportion)"), position = position_fill(0.7), vjust = 0.5, angle = 90, size = 3)+
      labs(y = "Proportion", x = column, fill = column2)

  }

  # Check if color_vector is provided and not NULL
  if (!is.null(color_vector)) {
    plot <- plot + scale_fill_manual(values = color_vector)
  }

  if (is.null(Legend_Title)) {
    plot <- plot  + labs(fill = "")
  }else{
    plot <- plot + labs(fill = Legend_Title)
  }

  if (is.null(Note_Title)) {
    plot <- plot # + labs(title = "")
  }else{
    plot <- plot + labs(title = paste(Note_Title, column2, ": Count based on", column))
  }

  plot <- plot +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          plot.title = element_text(size = 9, hjust = 0.5),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 12),
          aspect.ratio = aspect_ratio)

  return(plot)
}


# Plot_PvStat_pvalue_actual_cell_type <- plot_histogram(metadata, 'Actual Cell Type', 'PvStat', Note_Title = Note_TitleACT)
# Plot_Class_PvROC_clusters <- plot_histogram(metadata, 'seurat clusters', 'ClassificationROC', Note_Title = Note_TitleACT,color_vector = color_class)


################################################################################
# Helper function for plotting histogram with various settings
plot_histograms_for_vars <- function(data, var1, var2, Note_Title = "", color_vector = NULL) {
  list(
    default = plot_histogram(data, var1, var2, Note_Title = Note_Title, color_vector = color_vector, Legend_Title = ""),
    stack = plot_histogram(data, var1, var2, position_type = "stack", Note_Title = Note_Title, color_vector = color_vector, Legend_Title = ""),
    prop = plot_histogram(data, var1, var2, type = "proportion", Note_Title = Note_Title, color_vector = color_vector, Legend_Title = "")
  )
}

# Plot_PvStat_pvalue_actual_cell_type <- plot_histogram(metadata, 'Actual Cell Type', 'PvStat', Note_Title = Note_TitleACT)
# Plot_Class_PvROC_clusters <- plot_histogram(metadata, 'seurat clusters', 'ClassificationROC', Note_Title = Note_TitleACT,color_vector = color_class)


