# 整合各個子圖
combine_Figs <- function(plot_list, Set_ncol = 3) {
  plots_to_combine <- list()

  # 遍歷plot_list的每一項
  for (i in seq_along(plot_list)) {
    if (!is.null(plot_list[[i]])) {
      for (type in c( "stack", "prop", "default")) {
        # 檢查plot是否為NULL
        if (!is.null(plot_list[[i]][[type]])) {
          plots_to_combine[[paste0(i, "_", type)]] <- plot_list[[i]][[type]]
        }
      }
    }
  }

  # 確保至少有一個圖存在才進行組合
  if (length(plots_to_combine) > 0) {
    do.call(grid.arrange, c(plots_to_combine, ncol = Set_ncol))
  } else {
    message("No plots to combine!")
  }
}
