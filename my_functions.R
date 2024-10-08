

# Create effects plots and scatter points below

my_theme_square <- function() {
  theme_minimal(base_size = 8) +
    theme(
      #aspect.ratio = 1, 
      axis.ticks.y = element_line(),
      axis.ticks.x = element_line(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.background = element_rect(fill = "white", colour = "black"),
      legend.position = "bottom",
      axis.title.y = element_blank()
    ) 
}

# Create effect plot function with additional arguments to select columns from avg_data
create_effect_plot <- function(data, 
                               avg_data, 
                               x_col = "tmp_z_lag1", 
                               y_col = "sum_ips", 
                               line_color = "blue", 
                               x_title = "X-axis", 
                               y_title = "Y-axis", my_title = '',
                               x_annotate = 0, lab_annotate = "lab ann") {
  
  x_col <- ensym(x_col)
  y_col <- ensym(y_col)
  
  p <- ggplot() +
    geom_point(data = avg_data, aes(x = !!x_col, y = !!y_col), col = "gray60", alpha = 0.3) +
    geom_line(data = data, aes(x = x, y = predicted), color = line_color) +
    geom_ribbon(data = data, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = line_color) +
    labs(x = x_title,
         title = my_title,
         y = y_title) +
    my_theme_square() + 
    annotate("text", x = x_annotate, y = Inf, label = lab_annotate, hjust = 0.5, vjust = 1.5)
  
  return(p)
}

