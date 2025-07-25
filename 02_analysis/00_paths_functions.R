
# paths -------------------------------------
# define your path and folder
public_dir <- here("EU_forest_regeneration")

# Load cleaned input data

# final set of predictors, summary per plot
df_fin                 <- fread(file.path(public_dir, "data", "plot_level_predictors_clean.csv"))
# 
df_stem_species_class  <- fread(file.path(public_dir, "data", "plot_level_stem_density_species_by_class.csv"))

# simulated data aggregated on plot level
df_sim_indicators      <- fread(file.path(public_dir, "data", "df_sim_indicators.csv"))

# Species climate suitability :

# Read species lookup tables
lookup_species_acronyms      <- read.csv(file.path(public_dir, "data", "lookup_species_acronyms.csv"), sep = ';')
lookup_species_field_wessely <- read.csv(file.path(public_dir, "data", "lookup_species_field_wessely.csv"), sep = ';')

# read table with Wessely species - from species present/absent of plots locations
species_climate_suitability_raw <- fread(
  file.path(public_dir, "data", "species_presence_clim_change.csv")
)



### Get a color scheme per species -------------------------

# Reverse the color palette and map to the species in the desired order
n_colors  <- 10  # Number of species
my_colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))(n_colors)  # Generate colors

# Reverse the color order to start with dark green
reversed_colors <- rev(my_colors)

species_colors <- c(
  piab = "#006837",
  fasy = "#229C52",
  quro = "#74C364",
  pisy = "#B7E075",
  soau = "#E9F6A2",
  acps = "#FEEDA2",
  potr = "#FDBE6E",
  abal = "#F67B49",
  besp = "#DA362A",
  lade = "#A50026"
)


# Print the color assignments for confirmation
print(species_colors)

v_top10_species <- c("piab", "fasy", "quro", "pisy", "soau", "acps", "potr", "abal", "besp", "lade")

# Assign colors to each species in the order of `top_species_site_share$Species`
species_colors <- setNames(
  reversed_colors,
  c("piab", "fasy", "quro", "pisy", "soau", "acps", "potr", "abal", "besp", "lade")
)

# update species labels
species_labels <- c(
  "piab" = "Picea abies",
  "fasy" = "Fagus sylvatica",
  "quro" = "Quercus robur/petraea",
  "pisy" = "Pinus sylvestris",
  "soau" = "Sorbus aucuparia",
  "acps" = "Acer pseudoplatanus",
  "potr" = "Populus tremula",
  "abal" = "Abies alba",
  "besp" = "Betula sp.",
  "lade" = "Larix decidua"
)




# Variables ---------------------------------

# total number of plots
n_total_plots = length(unique(df_fin$plot)) # 849

reference_period <- 1980:2015
drought_period   <- 2018:2020
study_period     <- 2018:2023


# functions 

### Custom function to return median and IQR for stat_summary
median_iqr <- function(y) {
  median_val <- median(y, na.rm = TRUE)
  iqr_val <- IQR(y, na.rm = TRUE)
  ymin <- median_val - (iqr_val / 2)
  ymax <- median_val + (iqr_val / 2)
  return(c(y = median_val, ymin = ymin, ymax = ymax))
}




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
# SIMPLIFY PLOTTING TEST
# Function to generate predictions and plots
create_plot <- function(model, term, data, title, 
                        x_label = term, y_label = y_lab, 
                        line_color = "blue", fill_color = "blue", 
                        scatter = TRUE, x_limit = NULL, scatter_y = "stem_regeneration") {
  
  # Generate predictions
  predicted <- ggpredict(model, terms = term)
  
  # Create the base plot with scatter points first if required
  plot <- ggplot(predicted, aes(x = x, y = predicted)) +
    # Add scatter points behind the line and ribbon
    { if (scatter) geom_jitter(data = data, aes_string(x = term, y = scatter_y), 
                              color = "grey80", alpha = 0.5, size = 0.5) } +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = fill_color) +
    geom_line(color = line_color, size = 1) +
    labs(title = title, x = x_label, y = y_label) +
    theme_classic2() +
    theme(text = element_text(size = 8),
          panel.border = element_rect(color = "black", size = 0.7, fill = NA))
  
  # Set x-axis limits if provided
  if (!is.null(x_limit)) {
    plot <- plot + scale_x_continuous(limits = x_limit)
  }
  
  return(plot)
}
# Function to create interaction plots with scatter points in the background
create_interaction_plot <- function(model, terms, title, data, 
                                    x_label = terms[1], y_label = y_lab) {
  
  # Generate predictions for the interaction
  predicted_interaction <- ggpredict(model, terms = terms)
  
  # Create the plot with scatter points behind the line and ribbon
  plot <- ggplot(predicted_interaction, aes(x = x, y = predicted)) +
    # Scatter points in the background
    geom_jitter(data = data, aes_string(x = terms[1], y = "stem_regeneration"), 
               color = "grey80", alpha = 0.5, size = 0.5) +
    # Line and ribbon in the foreground
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
    geom_line(aes(color = group), size = 1) +
    # Labels and theme
    labs(title = title, x = x_label, y = y_label, color = "Group", fill = "Group") +
    theme_classic2() +
    theme(text = element_text(size = 8),
          legend.position = c(0, 1), legend.justification = c(0, 1),
          panel.border = element_rect(color = "black", size = 0.7, fill = NA))
  
  return(plot)
}



# no scatter plot ------------------------------------------

create_plot_no_scatter <- function(model, term, data, title, 
                        x_label = term, y_label = y_lab, 
                        line_color = "blue", fill_color = "blue", 
                        scatter = TRUE, x_limit = NULL, scatter_y = "stem_regeneration") {
  
  # Generate predictions
  predicted <- ggpredict(model, terms = term)
  
  # Create the base plot with scatter points first if required
  plot <- ggplot(predicted, aes(x = x, y = predicted)) +
    # Add scatter points behind the line and ribbon
   # { if (scatter) geom_jitter(data = data, aes_string(x = term, y = scatter_y), 
  #                             color = "grey80", alpha = 0.5, size = 0.5) } +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = fill_color) +
    geom_line(color = line_color, size = 1) +
    labs(title = title, x = x_label, y = y_label) +
    theme_classic2() +
    theme(text = element_text(size = 8),
          panel.border = element_rect(color = "black", size = 0.7, fill = NA))
  
  # Set x-axis limits if provided
  if (!is.null(x_limit)) {
    plot <- plot + scale_x_continuous(limits = x_limit)
  }
  
  return(plot)
}
# Function to create interaction plots with scatter points in the background
create_interaction_plot_no_scatter <- function(model, terms, title, data, 
                                    x_label = terms[1], y_label = y_lab) {
  
  # Generate predictions for the interaction
  predicted_interaction <- ggpredict(model, terms = terms)
  
  # Create the plot with scatter points behind the line and ribbon
  plot <- ggplot(predicted_interaction, aes(x = x, y = predicted)) +
    # Scatter points in the background
   # geom_jitter(data = data, aes_string(x = terms[1], y = "stem_regeneration"), 
    #            color = "grey80", alpha = 0.5, size = 0.5) +
    # Line and ribbon in the foreground
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
    geom_line(aes(color = group), size = 1) +
    # Labels and theme
    labs(title = title, x = x_label, y = y_label, color = "Group", fill = "Group") +
    theme_classic2() +
    theme(text = element_text(size = 8),
          legend.position = c(0, 1), legend.justification = c(0, 1),
          panel.border = element_rect(color = "black", size = 0.7, fill = NA))
  
  return(plot)
}


