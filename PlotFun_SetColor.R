update_color <- function(methods, color_method) {
  # Ensure color_method is a named vector
  color_method <- setNames(color_method, names(color_method))

  # Find methods without a color
  missing_methods <- setdiff(methods, names(color_method))

  # Assign a random color to each missing method
  set.seed(123) # For reproducibility
  colors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
  random_colors <- sample(colors, length(missing_methods))

  # Add the new color assignments to color_method
  color_method[missing_methods] <- random_colors

  return(color_method)
}

# ## Example usage:
# Assuming data_10x_SCINA$Method is a character vector of method names
# And color_Method is a named vector with names as method names and values as colors
# updated_color_method <- update_color(data_10x_SCINA$Method, color_Method)
