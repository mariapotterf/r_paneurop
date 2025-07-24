# 6. Delayed vs advanced sites: Wilcox plot-----------------------------------------------

public_dir <- here("outData", "public")

# Load cleaned input data
df_fin <- fread(file.path(public_dir, "data", "plot_level_predictors_clean.csv"))
df_stem_species_class <- fread(file.path(public_dir, "data", "plot_level_stem_density_species_by_class.csv"))


# Variables
# total number of plots
n_total_plots = length(unique(df_fin$plot)) # 849

source('00_my_functions.R')


# Combine predictor variables with the grouping variable
vars_to_plot <- c(predictor_vars_sub, "adv_delayed")


### test differences between plots:
# Step 1: Reshape the data to long format
df_long_narrow <- df_fin %>%
  na.omit() %>% 
  dplyr::select(all_of(vars_to_plot)) %>%
  gather(key = "Variable", value = "Value", -adv_delayed) %>% 
  droplevels() %>% 
  mutate(Variable = factor(Variable, levels = unique(Variable))) # preserve teh order of factors

# subset just few variables for a plot
df_long_narrow_sub <- df_fin %>%
  na.omit() %>% 
  dplyr::select( 
    prcp,
    tmp,
    tmp_z,
    disturbance_severity,
    distance_edge,
    clay_extract,
    drought_spei1,
    adv_delayed) %>% 
  #dplyr::select(all_of(variables_to_plot), adv_delayed) %>% 
  gather(key = "Variable", value = "Value", -adv_delayed) %>% 
  droplevels() %>% 
  mutate(Variable = factor(Variable, levels = unique(Variable))) # preserve teh order of factors

# Step 2: Calculate median and IQR for each variable by regeneration status
summary_stats_narrow_sub <- df_long_narrow_sub %>%
  na.omit() %>%
  group_by(adv_delayed, Variable) %>%
  summarise(
    Mean = round(mean(Value, na.rm = TRUE), 2),
    SD = round(sd(Value, na.rm = TRUE), 2),
    Median = round(median(Value, na.rm = TRUE), 2),
    IQR = round(IQR(Value, na.rm = TRUE), 2)
  ) %>%
  arrange(adv_delayed, Variable)

# list groups to pairwise comparison
comparisons_3grps <- list(c("Delayed", "Intermediate"), 
                          c("Delayed", "Advanced"), 
                          c("Intermediate", "Advanced"))



# Define a reusable theme function
my_theme <- function() {
  theme(
    legend.position = 'none',
    text = element_text(size = 8),         # Increase text size for readability
    axis.text = element_text(size = 8),    # Axis tick labels
    axis.title = element_text(size = 8),   # Axis titles
    strip.text = element_text(size = 8),   # Facet labels
    legend.text = element_text(size = 8),  # Legend text
    plot.title = element_text(size = 8),   # Plot title
    strip.background = element_blank(),    # Remove the box around facet names
    strip.placement = "outside",           # Move facet labels outside the plot area
    panel.border = element_rect(colour = "black",
                                fill = NA, 
                                linewidth = 0.5)  # Add a square border around the plots
  )
}


### Supplement Wilcox:  3  groups -------------------


# for all of predictors: 

# Define short labels
label_map <- c("Delayed" = "Del.", "Intermediate" = "Int.", "Advanced" = "Adv.")

# Function to generate one plot with custom p-value brackets
make_boxplot_with_pvals <- function(var_name, y_label, y_limits, bracket_heights, tick_height = 0.2, text_offset = 0.2) {
  
  # Filter data
  df_sub <- df_long_narrow_sub %>% filter(Variable == var_name)
  
  # Compute p-values
  pvals <- compare_means(Value ~ adv_delayed, data = df_sub,
                         method = "wilcox.test", comparisons = comparisons_3grps) %>%
    mutate(p_rounded = formatC(p, digits = 3, format = "f"))
  
  # Calculate dynamic tick height and text offset based on y-axis range
  y_range <- diff(range(y_limits))
  tick_height_scaled <- tick_height * y_range
  text_offset_scaled <- text_offset * y_range
  
  
  # Bracket data
  brackets <- data.frame(
    x = c(1, 1, 2),
    xend = c(2, 3, 3),
    y = bracket_heights,
    label = pvals$p_rounded,
    tick = tick_height,
    text_offset = text_offset
  )
  
  # Mean summary
  summary_data <- subset(summary_stats_narrow_sub, Variable == var_name)
  
  # Boxplot
  p <- ggboxplot(df_sub, x = "adv_delayed", y = "Value",
                 fill = "adv_delayed",
                 palette = c("#A50026", "#FDAE61", "#006837"),
                 ylab = y_label, xlab = "", outlier.shape = NA, size = 0.2) +
    scale_x_discrete(labels = label_map) +
    geom_point(data = summary_data,
               aes(x = adv_delayed, y = Mean),
               shape = 21, fill = "red", color = "red", size = 1.5, inherit.aes = FALSE) +
    coord_cartesian(ylim = y_limits) +
    my_theme() +
    
    # Bracket lines and ticks
    geom_segment(data = brackets, aes(x = x, xend = xend, y = y, yend = y), linewidth  = 0.3) +
    geom_segment(data = brackets, aes(x = x, xend = x, y = y, yend = y - tick), linewidth  = 0.3) +
    geom_segment(data = brackets, aes(x = xend, xend = xend, y = y, yend = y - tick), linewidth  = 0.3) +
    geom_text(data = brackets, aes(x = (x + xend) / 2, y = y + text_offset, label = label), linewidth  = 3)
  
  return(p)
}


p1 <- make_boxplot_with_pvals(
  "prcp", "Precipitation [mm]", 
  c(300, 1530), c(1270, 1370, 1490),
  tick_height = 20,       # ~1.7% of y-range
  text_offset = 35        # ~2.5% of y-range
)
p1
p2 <- make_boxplot_with_pvals(
  "tmp", "Temperature [°C]", 
  c(6.2, 14.2), c(12.6, 13.2, 14),
  tick_height = 0.2,       # ~2.5%
  text_offset = 0.25       # ~3%
)
p2

p3 <- make_boxplot_with_pvals(
  "drought_spei1", "SPEI-1 Drought Index [dim.]", 
  c(-1.15, -0.65), c(-0.8, -0.75, -0.7),
  tick_height = 0.009,     # ~3%
  text_offset = 0.02       # ~4%
)
p3
p4 <- make_boxplot_with_pvals(
  "clay_extract", "Clay [%]", 
  c(5, 55), c(42, 47, 53),
  tick_height = 1.2,       # ~2.7%
  text_offset = 1.7        # ~3.3%
)
p4
p5 <- make_boxplot_with_pvals(
  "disturbance_severity", "Disturbance severity [dim.]", 
  c(0, 1.5), c(1.1, 1.25, 1.4),
  tick_height = 0.03,      # ~2%
  text_offset = 0.04       # ~2.7%
)
p5
p6 <- make_boxplot_with_pvals(
  "distance_edge", "Distance to edge [m]", 
  c(0, 155), c(120, 133, 147),
  tick_height = 1.9,       # ~2.2%
  text_offset = 4          # ~2.5%
)
p6


## Wilcox supplement 
# Combine all plots into a single figure
wilcox_plot_supplement<- ggarrange(p1, p2, p3,p4,p5, p6,
                                   ncol = 3, nrow = 2, labels = c("[a]", "[b]","[c]","[d]", '[e]','[f]'), #'[g]','[h]'
                                   font.label = list(size = 8, face = "plain"))

# Save the plot as an SVG file
ggsave(
  filename = file.path(public_dir, "figs", "figS_boxplot_wilcox_multiple_vars.png"),
  plot = wilcox_plot_supplement,
  device = "png",
  width = 7,
  height = 5.5,
  dpi = 300,
  bg = "white"
)


### FINAL Wilcox: remove 'Intermediate' category -----

summary_stats_narrow_sub2 <- summary_stats_narrow_sub %>%
  dplyr::filter(adv_delayed != "Intermediate") %>%
  droplevels()

df_long_narrow_sub2 <- df_long_narrow_sub %>%
  dplyr::filter(adv_delayed != "Intermediate") %>%
  droplevels()



# Update the function to exclude "Intermediate" and use full labels
plot_variable_box_fin <- function(var_name, 
                                  y_label, 
                                  y_limits, 
                                  label_y_values, 
                                  tick_height = 0.02, 
                                  text_offset = 0.02, 
                                  palette = c("#F67B49", "#74C364")) {
  
  # Filter data and exclude Intermediate
  data_filtered <- df_long_narrow_sub2 %>%
    filter(Variable == var_name, adv_delayed %in% c("Delayed", "Advanced"))
  
  # Summary stats
  summary_data <- subset(summary_stats_narrow_sub2, 
                         Variable == var_name & adv_delayed %in% c("Delayed", "Advanced"))
  
  # Compute Wilcoxon p-value
  pval <- compare_means(Value ~ adv_delayed, 
                        data = data_filtered, 
                        method = "wilcox.test") %>%
    pull(p) %>%
    formatC(digits = 3, format = "f")
  
  # Prepare bracket layer
  bracket <- data.frame(
    x = 1,
    xend = 2,
    y = label_y_values,
    tick = tick_height,
    text_offset = text_offset,
    label = pval
  )
  
  # Create plot
  p <- ggboxplot(data_filtered,
                 x = "adv_delayed", 
                 y = "Value", 
                 fill = "adv_delayed", 
                 palette = palette,
                 ylab = y_label, 
                 xlab = "",
                 outlier.shape = NA,
                 size = 0.2) +
    
    # Mean dots
    geom_point(data = summary_data, 
               aes(x = adv_delayed, y = Mean), 
               shape = 21, fill = "black", color = "black", size = 1.5, inherit.aes = FALSE) +
    
    # Bracket line
    geom_segment(data = bracket, 
                 aes(x = x, xend = xend, y = y, yend = y), size = 0.3) +
    
    # Ticks
    geom_segment(data = bracket, 
                 aes(x = x, xend = x, y = y, yend = y - tick), size = 0.3) +
    geom_segment(data = bracket, 
                 aes(x = xend, xend = xend, y = y, yend = y - tick), size = 0.3) +
    
    # P-value text (no "p =")
    geom_text(data = bracket, 
              aes(x = (x + xend) / 2, y = y + text_offset, label = label), size = 3) +
    
    my_theme() +
    coord_cartesian(ylim = y_limits)
  
  return(p)
}


# Precipitation
p.prcp.fin <- plot_variable_box_fin(
  var_name = "prcp",
  y_label = "Precipitation [mm]",
  y_limits = c(300, 1500),
  label_y_values = c(1400),
  tick_height = 20,
  text_offset = 35
)

# SPEI-1
p.spei1.fin <- plot_variable_box_fin(
  var_name = "drought_spei1",
  y_label = "SPEI-1 [dim.]",
  y_limits = c(-1.15, -0.65),
  label_y_values = c(-0.7),
  tick_height = 0.009,
  text_offset = 0.02
)



# Clay
p.clay.fin <- plot_variable_box_fin(
  var_name = "clay_extract",
  y_label = "Clay [%]",
  y_limits = c(5,40),
  label_y_values = c(37),
  tick_height = 0.6,
  text_offset = 1
)




# Combine all plots into a single figure
wilcox_fin <- ggarrange(p.prcp.fin,
                        p.spei1.fin,p.clay.fin,
                        ncol = 3, 
                        nrow = 1, 
                        labels = c("[a]", "[b]","[c]"), font.label = list(size = 8, face = "plain"))

# Display the final merged plot
wilcox_fin


# Save the plot as an SVG file
ggsave(
  filename = file.path(public_dir, "figs", "fig_boxplot_wilcox_final_vars.png"),
  plot = wilcox_fin,
  device = "png",
  width = 6,
  height = 3,
  dpi = 300,
  bg = "white"
)


### use bayes statistics: more robust for uneven sample sizes ----------------------------
library("BayesFactor")
# Compute Bayesian ANOVA for each variable
bayesian_results <- df_long_narrow_sub2 %>%
  as.data.frame() %>% 
  mutate(adv_delayed = factor(adv_delayed)) %>% 
  
  group_by(Variable) %>%
  summarise(
    BF = extractBF(anovaBF(Value ~ adv_delayed, data = cur_data()))$bf[1]  # Extract Bayes Factor
  )# %>% 
# arrange(decs(BF))

bayesian_results <- bayesian_results %>%
  mutate(interpretation = case_when(
    BF < 1 ~ "Evidence for null",
    BF < 3 ~ "Anecdotal",
    BF < 10 ~ "Moderate",
    BF < 30 ~ "Strong",
    BF < 100 ~ "Very strong",
    BF >= 100 ~ "Decisive"
  ))


sjPlot::tab_df(bayesian_results,
               show.rownames = FALSE,
               file = file.path(public_dir, "tables", "table_bayesian_group_differences.doc"),
               digits = 1)

