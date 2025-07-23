# share analysis FigShare

gc()
# Libs --------------------------------------------------------------------------


library(data.table)
library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)
library(lubridate)
library(stringr)
library(ggpubr)
library(ggridges)
library(purrr)
library(mgcv)
library(gratia)
library(ggeffects)
#library(GGally) # for pairwise comparions of the variables 1:1, also calculates the correlation and sinificance

library(sf)     # read coordinates
library(scales)  # add % for the plots
library(RColorBrewer)
library(cowplot)

# Make a color gradient: shaded within species
library(colorspace)   # For gradient and color manipulation
library(here)
library(ggalluvial)

source('00_my_functions.R')

# Input data -------------------------------------------------------------
public_dir <- here("outData", "public")

# Load cleaned input data
df_fin <- fread(file.path(public_dir, "data", "plot_level_predictors_clean.csv"))
df_stem_species_class <- fread(file.path(public_dir, "data", "plot_level_stem_density_species_by_class.csv"))


# Variables
# total number of plots
n_total_plots = length(unique(df_fin$plot)) # 849

source('00_my_functions.R')

## Set a global theme -------------------------------------------------------------------

theme_set(
  theme_classic() + 
    theme(
      legend.position = 'bottom',
      text = element_text(size = 4),         # Set all text to size 8
      axis.text = element_text(size = 7),    # Axis tick labels
      axis.title = element_text(size = 7),   # Axis titles
      strip.text = element_text(size = 7),   # Facet labels
      legend.text = element_text(size =7),  # Legend text
      plot.title = element_text(size = 7)    # Plot title
    )
)



# Analysis ------------------------------------------------------------
# desription of current regeneration state
# structure 
# composition
# 
# What is the current state of the forest regeneration?
#  stucture? 
#    - stem density
#    - vertical structure
#  Composition? 
#    - richness
#    - species occurence


# Drivers: 
# effect of disturbance size:
#   - higher patch size, less regeneration
#   - higher severity, low regenerations

# environmental condistins;
# soil: 
#    - more stem density at better soil conditions (higher depth, nutrient content)


# how many plots do not have any stem density present? not evemn mature treees?
n_total <- length(unique(df_fin$plot))



# Categorize richness into bins and calculate proportions
richness_summary <- df_fin %>%
  mutate(
    richness_category = case_when(
      richness == 0 ~ "0",
      richness == 1 ~ "1",
      richness == 2 ~ "2",
      richness == 3 ~ "3",
      richness >= 4 ~ ">4",
      TRUE ~ "Other"  # Just in case there are unexpected values
    )
  ) %>%
  mutate(richness_category = factor(richness_category, levels = c("0","1","2","3",">4"))) %>% 
  group_by(richness_category) %>%
  summarise(plot_count = n_distinct(plot)) %>%
  mutate(proportion = plot_count / sum(plot_count) * 100)  # Calculate proportions

richness_summary

# Create a stacked bar plot
p_bar_richness_groups <- ggplot(richness_summary, 
                                aes(x = 1, 
                                    y = proportion, 
                                    fill = richness_category)) +
  geom_bar(stat = "identity", color ="black" 
  ) +
  geom_text(aes(label = paste0(round(proportion, 1), "%", '(', plot_count, ')')), 
            position = position_stack(vjust = 0.5), size = 2) +  # Add percentage labels
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    title = "",
    x = NULL,
    y = "Share of plots [%]",
    fill = "Species\nrichness"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text since it's a single bar
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, face = "plain"),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5)
  )

p_bar_richness_groups


### vertical classes share -----------------------------------
# get from vertical claas share

# Categorize counts of vertical claases into bins and calculate proportions
vert_str_summary <- df_fin %>%
  mutate(
    vert_category = case_when(
      n_vertical == 0 ~ "0",
      n_vertical == 1 ~ "1",
      n_vertical == 2 ~ "2",
      n_vertical == 3 ~ "3",
      TRUE ~ "Other"  # Just in case there are unexpected values
    )
  ) %>%
  mutate(vert_category = factor(vert_category, levels = c("0","1","2","3"))) %>% 
  group_by(vert_category) %>%
  summarise(plot_count = n_distinct(plot)) %>%
  mutate(proportion = plot_count / sum(plot_count) * 100)  # Calculate proportions

# Create a stacked bar plot
p_bar_vertical_groups <- ggplot(vert_str_summary, aes(x = 1, y = proportion, 
                                                      fill = vert_category)) +
  geom_bar(stat = "identity", color = "black"
  ) +
  geom_text(aes(label = paste0(round(proportion, 1), "%", '(', plot_count, ')')), 
            position = position_stack(vjust = 0.5), size = 2) +  # Add percentage labels
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    title = "",
    x = NULL,
    y = "Proportion (%)",
    fill = "Count vertical\nlayers"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text since it's a single bar
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, face = "plain"),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5)
  )


ggarrange(p_bar_vertical_groups, p_bar_richness_groups,  
          labels = c("[a]", "[b]"),  # Add labels
          ncol = 3,
          font.label = list(size = 10, face = "plain"))


## Species composition: Tables --------------

df_stem_sp_sum <- df_stem_species_class %>% 
  group_by(plot, Species) %>% 
  # summarize across vertical groups
  summarise(sum_stem_density = sum(stem_density, na.rm = t)) #%>% 


### Share of plots by species -----------------------------------------------------

# Count the unique plots (plots) for each species
species_plot_share <- df_stem_sp_sum %>%
  dplyr::filter(sum_stem_density >0) %>% 
  
  group_by(Species) %>%
  summarise(plot_count = n_distinct(plot)) %>%
  ungroup()


# Add a column for the share (percentage) of plots
species_plot_share <- species_plot_share %>%
  mutate(share_of_plots = (plot_count / n_total) * 100)

# Sort the results in descending order of share
species_plot_share <- species_plot_share %>%
  arrange(desc(share_of_plots))


top_species_plot_share <- species_plot_share %>% 
  arrange(desc(share_of_plots)) %>%
  slice_head(n = 10) # Select the top 10 species


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


v_top10_species <- top_species_plot_share$Species 
# #c("piab", "fasy", "quro", "pisy", "soau", "acps", "potr", "abal", "besp", "lade")

# Assign colors to each species in the order of `top_species_plot_share$Species`
species_colors <- setNames(
  reversed_colors,
  v_top10_species
)







## Create a horizontal bar plot -------------------------------------------------
p_plots_share_species <- species_plot_share %>%
  arrange(desc(share_of_plots)) %>%
  slice_head(n = 10) %>% # Select the top X rows
  mutate(Species = factor(Species, levels = unique(Species))) %>% 
  ggplot(aes(x = share_of_plots, y = reorder(Species, share_of_plots))) +
  geom_bar(stat = "identity", aes(fill = Species), #alpha = 0.8#, 
  ) + # Horizontal bars
  scale_fill_manual(values = species_colors) +  # Apply custom color palette
  labs(
    x = "\nShare of plots [%]",
    y = "",
    title = ""
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = '',
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title = element_text(size = 8),
    plot.title = element_text(size = 8, hjust = 0.5)
  )

p_plots_share_species



## Density plot of stem density  ----------------------------------

# Calculate median for each species and reorder the factor levels
df_stem_sp_sum_ordered <- df_stem_sp_sum %>%
  dplyr::filter(sum_stem_density >0) %>% 
  dplyr::filter(Species %in% v_top10_species ) %>%  #top_species_overall
  dplyr::group_by(Species) %>%
  dplyr::mutate(median_stem_density = median(sum_stem_density, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>%
  mutate(Species = factor(Species, levels = rev(v_top10_species))) # Set custom order
# dplyr::mutate(Species = reorder(Species, median_stem_density))  # Reorder species by median stem density

my_species_levels <-  levels(df_stem_sp_sum_ordered$Species)

# Add a log-transformed column for sum_stem_density
df_stem_sp_sum_ordered <- df_stem_sp_sum_ordered %>%
  mutate(log_sum_stem_density = log10(sum_stem_density + 1))  # Adding 1 to avoid log(0)

p_stem_density_species <- df_stem_sp_sum_ordered %>%
  ggplot(aes(x = log_sum_stem_density, y = Species, group = Species)) +
  geom_density_ridges(
    aes(fill = Species), 
    alpha = 1, 
    color = 'NA', 
    scale = 0.9 # Adjust the vertical scale of the ridges
  ) +
  scale_fill_manual(values = species_colors) +
  stat_summary(
    aes(x = log_sum_stem_density, fill = Species),  # Add fill aesthetic for inner color
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
    fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
    geom = "pointrange", 
    color = 'black'      ,  # Black outline for points
    shape = 21,  # Shape 21 is a circle with a fill and border
    size = 0.5,
    linewidth = 0.2,
    stroke = 0.2,
    position = position_nudge(y = 0.2)  # No vertical nudge for alignment
  ) +
  theme_classic() +
  labs(
    title = "",
    x = expression("\nStem density (log"[10]*") [n ha"^-1*"]"),
    # x = "\nStem density (log10) [#/ha]",
    y = ""
  ) +
  scale_x_continuous(
    labels = math_format(10^.x)  # Format x-axis labels as 10^3, 10^4, etc.
  ) +
  
  scale_y_discrete(expand = expansion(add = c(0.2, 0.2))) + # Adjust y-axis padding
  theme_classic(base_size = 8) +
  theme(
    axis.text = element_text(size = 8),    # Axis tick labels
    axis.title = element_text(size = 8),   # Axis titles
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    legend.position = 'none'               # Hide the legend
  )


# % of recorded stems of juveniles and saplings ---------------------
df_stem_species_class %>% 
  group_by(VegType) %>% 
  summarise(sum_vegType = sum(stem_density, na.rm = T)) %>% 
  mutate(total_stem_density = sum(sum_vegType),
         share = sum_vegType/total_stem_density)


## get % of stems per saplings/juveniles -----------------------------------------
# Summarize the total stem density per species 
species_composition_sapl_juv <- df_stem_species_class %>%
  dplyr::filter(Species %in% v_top10_species) %>% 
  dplyr::filter(VegType != "Mature") %>% 
  group_by(plot, Species,VegType) %>%
  summarize(sum_regeneration = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the stem density er regeneration class and share of each species
species_composition_sapl_juv <- 
  species_composition_sapl_juv %>%
  group_by(plot) %>%
  mutate(total_stems = sum(sum_regeneration),  # Total stem density in each climate class
         share = (sum_regeneration / total_stems) * 100) %>%  # Calculate percentage share
  ungroup() %>% 
  filter(!is.nan(share) & share > 0)

species_composition_sapl_juv_med <- species_composition_sapl_juv %>% 
  ungroup(.) %>% 
  group_by(Species, VegType) %>% 
  summarise(med_share = median(share),
            mean_share = mean(share))

# Transform data to make saplings negative for diverging plot
species_composition_plot_data <- species_composition_sapl_juv_med %>%
  mutate(
    med_share = ifelse(VegType == "Saplings", -med_share, med_share) # Saplings are negative
  ) %>% 
  mutate(Species = factor(Species, 
                          levels = rev(v_top10_species))) # Set custom order




# Calculate median and IQR for each species and VegType
species_composition_sapl_juv_summary <- species_composition_sapl_juv %>% 
  group_by(Species, VegType) %>% 
  summarise(
    med_share = median(share, na.rm = TRUE),
    iqr_low = quantile(share, 0.25, na.rm = TRUE),
    iqr_high = quantile(share, 0.75, na.rm = TRUE),
    iqr = IQR(share, 0.75, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    med_share = ifelse(VegType == "Saplings", -med_share, med_share),  # Saplings are negative
    iqr_low = ifelse(VegType == "Saplings", -iqr_low, iqr_low),       # Adjust IQR for saplings
    iqr_high = ifelse(VegType == "Saplings", -iqr_high, iqr_high),    # Adjust IQR for saplings
    Species = factor(Species, levels = rev(v_top10_species)) # Custom species order
  )

species_composition_sapl_juv_summary %>% 
  group_by(VegType) %>% 
  summarize(median_iqr = median(iqr),
            min_iqr = min(iqr),
            max_iqr = max(iqr))



# Calculate overall median share for Juveniles and Saplings
overall_med_share <- species_composition_sapl_juv_summary %>%
  group_by(VegType) %>%  # Group by Vegetation Type (Juveniles/Saplings)
  summarize(
    overall_med_share = median(med_share, na.rm = TRUE),
    IRQ_med_share = IQR(med_share, na.rm = TRUE),
    min_med_share = min(med_share, na.rm = TRUE),
    max_med_share = max(med_share, na.rm = TRUE)# Calculate median
  )



## adjust colors -----------------------------------------------------------------

# Generate separate colors for Juveniles (darker) and Saplings (lighter)
species_colors_juveniles <- lapply(species_colors, function(color) colorspace::darken(color, 
                                                                                      amount = 0.2)) %>% unlist()

species_colors_saplings <- lapply(species_colors, function(color) colorspace::lighten(color, 
                                                                                      amount = 0.2)) %>% unlist()

# Create a named vector combining Juveniles and Saplings
species_colors_combined <- c(
  setNames(species_colors_juveniles, paste0(names(species_colors), "_Juveniles")),
  setNames(species_colors_saplings, paste0(names(species_colors), "_Saplings"))
)

# Add a combined key to species_composition_sapl_juv_summary for fill mapping
species_composition_sapl_juv_summary <- species_composition_sapl_juv_summary %>%
  mutate(fill_key = paste0(Species, "_", VegType))

# Diverging bar chart with species-specific colors
p_share_vertical_species <- 
  ggplot(species_composition_sapl_juv_summary, aes(x = med_share, y = Species, fill = fill_key)) +
  geom_bar(stat = "identity", position = "identity", alpha = 1) + # Horizontal bars
  geom_errorbarh(aes(xmin = iqr_low, xmax = iqr_high), 
                 height = 0.1, 
                 color = "black",                   ,
                 linewidth = 0.2) + # Add IQR error bars
  scale_x_continuous(
    labels = abs, # Show positive labels
    name = "Stem density\nMedian share [%]",
  ) +
  scale_fill_manual(values = species_colors_combined) + # Use combined color palette
  labs(
    x = "Stem density\nMedian share [%]",
    y = "",
    title = ""
  ) +
  geom_vline(xintercept = 0, linetype = "solid", 
             color = "black",
             linewidth = 0.2) + #
  theme_classic(base_size = 8) +
  annotate("text", x = -60, 
           #y = length(levels(species_composition_sapl_juv_summary$Species)) + 0.5,
           y = 1,
           label = "Saplings", color = "gray30", size = 3, fontface = "plain") + # Add Saplings label
  annotate("text", x = 50, 
           #y = length(levels(species_composition_sapl_juv_summary$Species)) + 0.5,
           y = 1,
           label = "Juveniles", color = "gray30", size = 3, fontface = "plain") + # Add Juveniles label
  
  theme(
    axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 8),
    legend.position = "none", # Remove legend as it will clutter
    plot.title = element_text(size = 8, hjust = 0.5)
  )



# Remove x- and y-axis labels from individual plots
p1 <- p_plots_share_species +
  scale_y_discrete(labels = species_labels) + # Replace y-axis labels with full Latin names
  
  theme(
    axis.text.y = element_text(face = "italic", size = 8)
  )
p2 <- p_stem_density_species + theme(axis.ticks.y = element_blank(),
                                     axis.text.y = element_blank(), 
                                     axis.title.y = element_blank())
p3 <- p_share_vertical_species + theme(axis.ticks.y = element_blank(), 
                                       axis.text.y = element_blank(), 
                                       axis.title.y = element_blank())


# Combine plots with aligned axes
combined_plot <- plot_grid(
  p1,p2,p3,
  align = "h",  # Align vertically
  axis = "b",
  ncol = 3 ,     # Three plots in a row
  rel_widths = c(1.4, 1, 1), # Adjust relative widths if necessary
  labels = c("[a]", "[b]", "[c]"), # Add labels
  label_size = 8,           # Adjust label size
  label_fontface = "plain"   # Use plain font for labels
)

# Save the combined plot as an image
ggsave(
  filename = file.path(public_dir, "figs", "Fig1.png"),    # File name (change extension for different formats, e.g., .pdf)
  plot = combined_plot,             # Plot object to save
  width = 7,                       # Width of the saved plot in inches
  height = 3.5,                       # Height of the saved plot in inches
  dpi = 300,                        # Resolution (dots per inch)
  units = "in"                      # Units for width and height
)


## Bar plot with error plot ---------------------------------------------------------

# Add a log-transformed column for sum_stem_density
df_stem_species_class <- df_stem_species_class %>%
  mutate(log_stem_density = log10(stem_density + 1))  # Adding 1 to avoid log(0)


# simplify naming and make new categories for proper stem density
df_stem_species_class <- df_stem_species_class %>%
  mutate(VegType_acc = case_when(
    VegType == "Mature" ~ "mat",
    VegType == "Juveniles" ~ "juv",
    VegType == "Saplings" ~ "sap",
    TRUE ~ as.character(VegType)  # Keep any other VegType values as they are
  )) %>% 
  mutate(species_VegType = paste(Species, VegType_acc, sep = '_')) %>% 
  mutate(VegType_acc = factor(VegType_acc, 
                              levels = c('mat', 'juv', 'sap'))) 




## Tableplot:  Stem density per species and vertical class (): --------------------------------

species_lookup <- c(
  abal = "Abies alba",
  acca = "Acer campestre",
  acpl = "Acer platanoides",
  acps = "Acer pseudoplatanus",
  algl = "Alnus glutinosa",
  alin = "Alnus incana",
  besp = "Betula sp.",
  cabe = "Carpinus betulus",
  casa = "Castanea sativa",
  fasy = "Fagus sylvatica",
  frex = "Fraxinus excelsior",
  fror = "Fraxinus ornus",
  ilaq = "Ilex aquifolium",
  juni = "Juniperus sp.",
  jure = "Juglans regia",
  lade = "Larix decidua",
  osca = "Ostrya carpinifolia",
  piab = "Picea abies",
  pist = "Pinus strobus",
  pisy = "Pinus sylvestris",
  posp = "Populus sp.",
  potr = "Populus tremula",
  prav = "Prunus avium",
  psme = "Pseudotsuga menziesii",
  quro = "Quercus robur/rubra",
  qusp = "Quercus sp.",
  rops = "Robinia pseudoacacia",
  saca = "Salix caprea",
  sasp = "Salix sp.",
  soar = "Sorbus aria",
  soau = "Sorbus aucuparia",
  soto = "Sorbus torminalis",
  taba = "Tilia cordata",
  tisp = "Tilia sp.",
  ulsp = "Ulmus sp."
)

# Get a summary table if species is present:
summary_stem_dens_VegType <- df_stem_species_class %>%
  #  dplyr::filter(Species %in% v_top10_species) %>%
  dplyr::filter(stem_density > 0) %>%
  #  mutate(Species = factor(Species, levels = v_top10_species)) %>%
  group_by(Species, VegType_acc) %>%
  summarise(min = min(stem_density, na.rm =T),
            max = max(stem_density, na.rm =T),
            mean     = mean(stem_density, na.rm =T),
            sd = sd(stem_density, na.rm =T),
            Median = median(stem_density, na.rm =T),
            IQR = IQR(stem_density)#,
            #Q1 = quantile(stem_density, 0.25, na.rm =T),  # 25th percentile
            # Q3 = quantile(stem_density, 0.75, na.rm =T)   # 75th percentile
  ) %>%
  arrange(Species, VegType_acc)

summary_stem_dens_VegType

# # Get a summary table:
summary_stem_dens_spec <- df_stem_species_class %>%
  #  dplyr::filter(Species %in% v_top10_species) %>%
  dplyr::filter(stem_density > 0) %>%
  # mutate(Species = factor(Species, levels = v_top10_species)) %>%
  group_by(Species) %>%
  summarise(
    min = min(stem_density, na.rm =T),
    max = max(stem_density, na.rm =T),
    mean     = mean(stem_density, na.rm =T),
    sd = sd(stem_density, na.rm =T),
    Median = median(stem_density, na.rm =T),
    IQR = IQR(stem_density)#,
    
  ) %>%
  mutate(Species = recode(Species, !!!species_lookup)) %>% 
  arrange(Species)



# Display the summary table
summary_stem_dens_spec

# print ourt

# Reorder and concatenate columns in the desired format
summary_stem_dens_spec_formated <- summary_stem_dens_spec %>%
  mutate(
    #min_max = paste(min, max, sep = " - "),
    mean_sd = paste0(round(mean, 0), " ± ", round(sd, 0)),
    median_iqr = paste0(Median, " (", IQR, ")")
  ) %>%
  dplyr::select(Species, min, max, mean_sd, median_iqr)

# Display the updated table
summary_stem_dens_spec_formated

# Export the table using sjPlot
sjPlot::tab_df(
  summary_stem_dens_spec_formated, 
  file = file.path(public_dir, "tables", "table_species_stem_density_summary.doc"),
  title = "Summary of Stem Density per Species",
  show.rownames = FALSE
)



## Species composition: vertical class  ------------------------------------
# Summarize the data by VegType and Species
mean_stems_vertical <- df_stem_species_class %>%
  group_by(VegType, Species) %>%
  summarize(mean_stems = mean(stem_density, na.rm = TRUE)) %>%
  ungroup()

# Calculate total stems per species and reorder Species factor
species_order <- mean_stems_vertical %>%
  group_by(Species) %>%
  summarize(total_stems = sum(mean_stems)) %>%
  arrange(total_stems)

# Reorder Species based on total stem density (ascending, so most prevalent species will be in the back)
mean_stems_vertical <- mean_stems_vertical %>%
  mutate(Species = factor(Species, levels = species_order$Species))

# Reorder VegType factor to have "Saplings", "Juveniles", "Mature"
mean_stems_vertical <- mean_stems_vertical %>%
  mutate(VegType = factor(VegType, levels = c("Saplings", "Juveniles", "Mature")))



#### Species composition : Richness ----------------------
# Summarize the total stem density per species
species_composition <- df_stem_species_class %>%
  group_by(Species) %>%
  summarize(sum_stems = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the total stem density per climate class and the share of each species
species_composition <- species_composition %>%
  #  group_by(clim_class) %>%
  mutate(total_stems = sum(sum_stems),  # Total stem density in each climate class
         share = (sum_stems / total_stems) * 100) %>%  # Calculate percentage share
  ungroup()



## Structure ---------------------------------------------------------

### Vertical structure  --------------------------------------

# Group by plot and VegType, check if stem_density > 0
vert_class_presence_absence <- df_stem_species_class %>%
  group_by(plot, VegType_acc) %>%
  summarise(has_stem_density = sum(stem_density > 0, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  # Filter only where there is stem density present
  dplyr::filter(has_stem_density) %>%
  dplyr::select(plot, VegType_acc) %>% 
  mutate(Presence = 1) %>%  # Add a column with 1 indicating presence
  pivot_wider(names_from = VegType_acc, values_from = Presence, values_fill = 0)  # Convert to wide format, filling NAs with 0



# Now find plots where no VegType has stem_density > 0
# Find plots where no vegType has stem_density
all_plots <- df_stem_species_class %>%
  ungroup() %>% 
  dplyr::select(plot, VegType_acc) %>% 
  distinct(plot)

# Combine plots with and without stem density
vert_class_presence_absence_fin <- vert_class_presence_absence  %>% 
  right_join(all_plots)




# 5. Drivers ---------------------------------

# select main variables as predictors 
predictor_vars_sub <- c(
  
  # over 2018-2023
  "spei1",
  "spei12",
  "tmp", 
  "prcp", 
  "tmp_z", 
  "prcp_z", 
  
  # during 2018-2020
  "drought_tmp",
  "drought_prcp",
  "drought_spei1",
  "drought_spei12",
  
  # disturbance chars
  "distance_edge", 
  "disturbance_severity", # from RS 
 
  # plot info
  "sand_extract",
  "clay_extract", 
  "depth_extract", 
  "av.nitro",
 
  # seasonality: CV - over year
  "cv_t2m",
  "cv_tp" #,
)


# Subset the data 

df_stem_regeneration <- df_fin %>% 
  dplyr::select(all_of(c("plot", 
                         "stem_regeneration",
                         predictor_vars_sub,
                         "country_pooled", 
                         "x", "y"))) %>% 
  mutate(country_pooled = factor(country_pooled))


summary(df_stem_regeneration)


## Pre-diagnostics --------------------------------
### Univariate models to find teh best predictors  -----------------------------

# Define response and predictor variables
response_var <- "stem_regeneration" # Replace with your actual response variable name
predictor_vars <- predictor_vars_sub

# Create an empty data frame to store results
AIC_results_univariate <- data.frame(
  Predictor = character(),
  AIC = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each predictor variable
for (predictor in predictor_vars) {
  # Create the formula
  formula <- as.formula(paste(response_var, "~ s(", predictor, ", k = 3)"))
  
  # Fit the GAM model with Tweedie distribution
  model <- gam(formula, family = tw(), data = df_stem_regeneration )
  
  # Extract AIC
  aic <- AIC(model)
  
  # Store the results
  AIC_results_univariate <- rbind(AIC_results_univariate, data.frame(Predictor = predictor, AIC = aic))
}

# Sort results by AIC (lower is better)
AIC_results_univariate <- AIC_results_univariate[order(AIC_results_univariate$AIC), ]

# Display the results
View(AIC_results_univariate)

sjPlot::tab_df(AIC_results_univariate,
               show.rownames = FALSE,
               file = file.path(public_dir, "tables", "table_univariate_GAM_AIC_comparison.doc"),
               digits = 1) 


### Spearman correlation coefficient ---------------------------------------------

# check correlation between predictors
library(corrplot)

# Select the relevant predictors from your data frame
predictors <- df_stem_regeneration %>%
  # dplyr::select(where(is.numeric))
  dplyr::select(all_of(predictor_vars_sub))

# Calculate the correlation matrix
correlation_matrix <- cor(predictors, use = "complete.obs", method = "spearman")

# Rename variables in the correlation matrix
correlation_matrix_renamed <- correlation_matrix

# Rename rows
rownames(correlation_matrix_renamed) <- rownames(correlation_matrix_renamed) %>%
  stringr::str_replace("^cv_t2m$", "cv_tmp") %>%
  stringr::str_replace("^cv_tp$", "cv_prcp")

# Rename columns
colnames(correlation_matrix_renamed) <- colnames(correlation_matrix_renamed) %>%
  stringr::str_replace("^cv_t2m$", "cv_tmp") %>%
  stringr::str_replace("^cv_tp$", "cv_prcp") %>% 
  stringr::str_replace("_extract$", "_content")


# Export with high resolution and clearer layout
png(file.path(public_dir, "figs", "fig_correlation_matrix_predictors.png"), 
    width = 2000, height = 2000, res = 300)

corrplot(correlation_matrix_renamed, 
         method = "color",
         type = "upper",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.col = "black",         # Label color
         tl.cex = 0.8,             # Smaller label font
         tl.srt = 45,              # Diagonal tilt (more readable)
         #addCoef.col = NA,         # Omit numbers for clarity
         # number.cex = 0.6,         # If you later add coefficients
         mar = c(0, 0, 1, 0)       # Reduce margins
)

dev.off()

## GAM(M)s------------

m_rnd <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(prcp,tmp, k = 5 ) +
    s(country_pooled, bs = "re") +
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration
)

m_rnd_te <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(country_pooled, bs = "re") +
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration
)
m_rnd_ti <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(country_pooled, bs = "re") +
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration
)

m_rnd_ti_fixed <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
       s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration
)

AIC(m_rnd_ti_fixed, m_rnd_ti)

### Diagnostics  ----------------

mgcv::summary.gam(m_rnd_ti_fixed)
gratia::appraise(m_rnd_ti_fixed)

#### Test for spatial autocorrelation:  

# Extract model residuals
m <- fin.m.reg.density
model_residuals <- residuals(m , type = "pearson")  

# Load necessary libraries
library(spdep)

# Create coordinates matrix
coords <- cbind(m$x, m$y)

# Create a spatial neighbors object (e.g., using k-nearest neighbors)
# Adjust k based on the density and distribution of your data points
nb <- knn2nb(knearneigh(coords, k = 5))

# Convert neighbors list to a weights list
listw <- nb2listw(nb, style = "W")
# Perform Moran's I test
moran_test <- moran.test(model_residuals, listw)
print(moran_test)


### Save model --------------------------------------------------
fin.m.reg.density <- m_rnd_ti_fixed

save(fin.m.reg.density, df_stem_regeneration,
     file = file.path(public_dir, "model", "model_stem_regeneration_density.RData"))



#### Relative contribution of factors  ------------------------------------
### drop one variable at time:  

# Fit full model
m_full <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5) +
    ti(prcp, tmp, k = 5) +
    s(x, y),
  family = tw(),
  method = "REML", select = TRUE,
  data = df_stem_regeneration
)
# Store logLik and AIC
loglik_full <- logLik(m_full)
aic_full <- AIC(m_full)

# Build and compare models one by one
m_drop_prcp <- gam(stem_regeneration ~ 
                     s(tmp, k = 5) + s(distance_edge, k = 5) +
                     s(disturbance_severity, k = 5) + s(clay_extract, k = 5) +
                     s(av.nitro, k = 5) + ti(disturbance_severity, distance_edge, k = 5) +
                     ti(prcp, tmp, k = 5) + s(x, y),
                   family = tw(), method = "REML",  select = TRUE, data = df_stem_regeneration)

m_drop_tmp <- gam(stem_regeneration ~ 
                    s(prcp, k = 5) + s(distance_edge, k = 5) +
                    s(disturbance_severity, k = 5) + s(clay_extract, k = 5) +
                    s(av.nitro, k = 5) + ti(disturbance_severity, distance_edge, k = 5) +
                    ti(prcp, tmp, k = 5) + s(x, y),
                  family = tw(), method = "REML",  select = TRUE,data = df_stem_regeneration)

m_drop_dist_edge <- gam(stem_regeneration ~ 
                          s(prcp, k = 5) + s(tmp, k = 5) +
                          s(disturbance_severity, k = 5) + s(clay_extract, k = 5) +
                          s(av.nitro, k = 5) + ti(disturbance_severity, distance_edge, k = 5) +
                          ti(prcp, tmp, k = 5) + s(x, y),
                        family = tw(), method = "REML",  select = TRUE, data = df_stem_regeneration)

m_drop_dist_sev <- gam(stem_regeneration ~ 
                         s(prcp, k = 5) + s(tmp, k = 5) +
                         s(distance_edge, k = 5) + s(clay_extract, k = 5) +
                         s(av.nitro, k = 5) + ti(prcp, tmp, k = 5) + s(x, y),
                       family = tw(), method = "REML",  select = TRUE, data = df_stem_regeneration)

m_drop_prcp_tmp <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5) +
    # ti(prcp, tmp, k = 5) +
    s(x, y),
  family = tw(),
  method = "REML",
  select = TRUE,
  data = df_stem_regeneration
)

# Compare models
results <- data.frame(
  model = c("m_rnd_ti_fixed", "drop_prcp", "drop_tmp", "drop_dist_edge", "drop_dist_sev", 'm_drop_prcp_tmp'),
  AIC = c(AIC(m_rnd_ti_fixed), AIC(m_drop_prcp), AIC(m_drop_tmp), AIC(m_drop_dist_edge), AIC(m_drop_dist_sev), AIC(m_drop_prcp_tmp)),
  logLik = c(logLik(m_rnd_ti_fixed), logLik(m_drop_prcp), logLik(m_drop_tmp), logLik(m_drop_dist_edge), logLik(m_drop_dist_sev), logLik(m_drop_prcp_tmp)),
  R2_adj = c(summary(m_rnd_ti_fixed)$r.sq, summary(m_drop_prcp)$r.sq, summary(m_drop_tmp)$r.sq,
             summary(m_drop_dist_edge)$r.sq, summary(m_drop_dist_sev)$r.sq, summary(m_drop_prcp_tmp)$r.sq),
  deviance_expl = c(summary(m_rnd_ti_fixed)$dev.expl, summary(m_drop_prcp)$dev.expl, summary(m_drop_tmp)$dev.expl,
                    summary(m_drop_dist_edge)$dev.expl, summary(m_drop_dist_sev)$dev.expl, summary(m_drop_prcp_tmp)$dev.expl)
)

results
# 
# > results
# model      AIC    logLik     R2_adj deviance_expl
# 1  m_rnd_ti_fixed 16069.08 -8013.875 0.09310789    0.10887421
# 2       drop_prcp 16072.64 -8011.886 0.09795147    0.11237287
# 3        drop_tmp 16070.49 -8013.018 0.09464378    0.11037273
# 4  drop_dist_edge 16073.63 -8015.650 0.09262081    0.10577277
# 5   drop_dist_sev 16076.56 -8019.264 0.09183355    0.09942951
# 6 m_drop_prcp_tmp 16070.75 -8015.014 0.09329267    0.10687028

### Low vs high severity plots ----------------------------

# Calculate the percentage of plots with >90% severity
prop_high_severity <- mean(df_stem_regeneration$disturbance_severity > 0.9, na.rm = TRUE) * 100
prop_high_severity

# Subset: Only plots with disturbance severity < 90%
df_low_severity <- df_stem_regeneration %>%
  dplyr::filter(disturbance_severity < 0.7)

nrow(df_low_severity)

# Refit the same GAM on this subset
m_lowsev <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5) +
    ti(prcp, tmp, k = 5) +
    s(x, y),
  family = tw(),
  method = "REML",
  select = TRUE,
  data = df_low_severity
)
plot(m_lowsev, pages = 1)
plot(m_rnd_ti_fixed, pages = 1)


summary(m_lowsev)
summary(m_rnd_ti_fixed)


#low vs full severity plots

# 1. Predict from the full model (across full severity range)
pred_full <- ggpredict(m_rnd_ti_fixed, terms = "disturbance_severity [0:1 by=0.01]") # 
pred_full$group <- "Full model (n = 849)"

# 2. Predict from the low-severity model (limited to 0–0.7 range)
pred_low <- ggpredict(m_lowsev, terms = "disturbance_severity [0:0.7 by=0.01]") #  
pred_low$group <- "Low-severity subset (n = 133)"

# 3. Combine predictions
pred_combined <- bind_rows(as.data.frame(pred_full), as.data.frame(pred_low)) %>%
  mutate(disturbance_severity_pct = x * 100)

# 4. Plot both prediction lines + confidence bands + raw data
p_low_high_severity <- 
  ggplot(pred_combined, aes(x = disturbance_severity_pct, y = predicted/1000, 
                            color = group,
                            linetype = group)) +
  geom_ribbon(aes(ymin = conf.low/1000, 
                  ymax = conf.high/1000, fill = group), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Full model (n = 849)" = "steelblue",
                                "Low-severity subset (n = 133)" = "darkorange"),
                     guide = guide_legend(nrow = 2)) +
  scale_fill_manual(values = c("Full model (n = 849)" = "steelblue",
                               "Low-severity subset (n = 133)" = "darkorange"),
                    guide = guide_legend(nrow = 2)) +
  scale_linetype_manual(values = c("Full model (n = 849)" = "solid",
                                   "Low-severity subset (n = 133)" = "dashed"),
                        guide = guide_legend(nrow = 2)) +
  
  labs(
    x = "Disturbance severity [%]",
    y = expression("Stem regeneration density [*1000 n ha"^-1*"]"),
    color = "", fill = "", linetype = "",
  ) +
  theme_classic2(base_size = 8) +
  theme(legend.position = "bottom")


p_low_high_severity

# Save the last plot
ggsave(p_low_high_severity,
       filename = file.path(public_dir, "figs", "figS_stem_density_by_disturbance_severity_subset_vs_full.png"),
       width = 4,
       height = 4,
       units = "in",
       dpi = 300  # High resolution for publication
)






### Export final drivers model: only fixed effectd  -----------------------------------------------------------------

# Identify random effects using the model's "smooth" component
smooth_terms <- summary(fin.m.reg.density)$s.table
# Identify random effects using the model's "smooth" component
smooth_terms_lowsev <- summary(m_lowsev)$s.table


# Extract the smooth terms labels and check which ones are random effects
random_effects_labels <- rownames(smooth_terms)[str_detect(rownames(smooth_terms), "country_pooled|clim_grid")]
random_effects_labels_lowsev <- rownames(smooth_terms_lowsev)[str_detect(rownames(smooth_terms_lowsev), "country_pooled|clim_grid")]


# Create a function to automatically label the terms in the summary output
create_labels <- function(term) {
  if (term %in% random_effects_labels) {
    return(paste(term, "(Random Effect)"))
  } else {
    return(term)
  }
}

# Apply the labeling function
pred_labels <- sapply(rownames(smooth_terms), create_labels)
pred_labels_lowsev <- sapply(rownames(smooth_terms_lowsev), create_labels)

# Save tab_model for full model
sjPlot::tab_model(
  fin.m.reg.density,
  show.re.var = TRUE,
  pred.labels = c("Intercept", pred_labels),
  dv.labels = paste0("Explained Deviance: ", round(100 * summary(fin.m.reg.density)$dev.expl, 2), "%"),
  file = file.path(public_dir, "tables", "table_model_stem_regeneration_density_full.doc")
)

# Save tab_model for low-severity subset model
sjPlot::tab_model(
  m_lowsev,
  show.re.var = TRUE,
  pred.labels = c("Intercept", pred_labels_lowsev),
  dv.labels = paste0("Explained Deviance: ", round(100 * summary(m_lowsev)$dev.expl, 2), "%"),
  file = file.path(public_dir, "tables", "table_model_stem_regeneration_density_low_severity.doc")
)

#### with random effects 

# Identify random effects using the model's "smooth" component
smooth_terms <- summary(m_rnd_ti)$s.table

# Extract the smooth terms labels and check which ones are random effects
random_effects_labels <- rownames(smooth_terms)[str_detect(rownames(smooth_terms), "country_pooled|clim_grid")]

# Create a function to automatically label the terms in the summary output
create_labels <- function(term) {
  if (term %in% random_effects_labels) {
    return(paste(term, "(Random Effect)"))
  } else {
    return(term)
  }
}

# Apply the labeling function
pred_labels <- sapply(rownames(smooth_terms), create_labels)

# Display the tab_model with automatic labeling
sjPlot::tab_model(
  m_rnd_ti,
  show.re.var = TRUE,
  pred.labels = c("Intercept", pred_labels),
  dv.labels = paste0("Explained Deviance: ", round(100 * summary(m_rnd_ti)$dev.expl, 2), "%"),
  file = file.path(public_dir, "tables", "table_model_stem_regeneration_density_full_with_random_effects.doc")
)




### Plot: Drivers  ---------------------------------------------------------------------------

y_lab = expression("Stem density [1000 n ha"^{-1}*"]")
y_lab_reg = expression("Reg. stem density [1000 n ha"^{-1}*"]")

m <- fin.m.reg.density #m_int_sev_edge_full_te_comb # m_int_res_edge_full_te_comb #  m_int_res_edge_full_te

# Generate predictions using ggpredict
summary(m)

# Interaction 1: Precipitation and Temperature
pred_prcp       <- ggpredict(m, terms = c("prcp"))
pred_tmp        <- ggpredict(m, terms = c("tmp"))
pred_tmp_prcp   <- ggpredict(m, terms = c("prcp", "tmp [8,10]"))
pred_dist_edge  <- ggpredict(m, terms = c("distance_edge[50:250]"))
pred_dist_sever <- ggpredict(m, terms = c("disturbance_severity"))
pred_clay       <- ggpredict(m, terms = c("clay_extract"))

# Convert all ggpredict objects to data.frames
pred_prcp        <- as.data.frame(pred_prcp)
pred_tmp         <- as.data.frame(pred_tmp)
pred_tmp_prcp    <- as.data.frame(pred_tmp_prcp)
pred_dist_edge   <- as.data.frame(pred_dist_edge)
pred_dist_sever  <- as.data.frame(pred_dist_sever)
pred_clay        <- as.data.frame(pred_clay)


my_theme_drivers <- theme(
  axis.title = element_text(size = 8),
  plot.title = element_text(hjust = 0.5, size = 8),       # Title size
  #axis.title.y = element_blank(),                   # Axis title size
  axis.text = element_text(size = 8),                    # Axis text size
  legend.key.size = unit(0.5, "cm"),                     # Legend key size
  legend.text = element_text(size = 8),                  # Legend text size
  legend.title = element_text(size = 8),                 # Legend title size
  legend.position = c(0.05, 0.9),
  legend.justification = c(0.05, 0.9)
)

#my_colors_interaction <- c("grey90", "red") 
my_colors_interaction <- c("#FDAE61", "#A50026") 
my_color_main_effects <- "grey90" # "#006837"    

# Plot precipitation
p.prcp <- ggplot(pred_prcp, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), alpha = 0.4,fill = my_color_main_effects) +
  labs(x = "Precipitation [mm]", 
       y = y_lab_reg, 
       title = "p<0.001") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p.prcp

# Plot precipitation
p.tmp <- ggplot(pred_tmp, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), alpha = 0.4,fill = my_color_main_effects) +
  labs(x = "Temperature [°C]",  
       y = y_lab_reg, 
       title = "p=0.004") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p.tmp



# Plot the first interaction
p1 <- 
  ggplot(pred_tmp_prcp, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, aes(color = group) ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000, fill = group), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  scale_fill_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  # theme_classic() +
  labs(x = "Precipitation [mm]", 
       y =y_lab_reg,  
       title = "p=0.006", 
       # linetype =  "Temperature [°C]"
  ) +
  my_theme_drivers

p1

# Plot distance to edge
p2 <- ggplot(pred_dist_edge, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), 
              alpha = 0.3,fill = my_color_main_effects) +
  labs(x = "Distance to edge [m]",  
       y = y_lab_reg, 
       title = "p=0.023") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p2


# 
p3 <- ggplot(pred_dist_sever, aes(x = x*100, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), fill =my_color_main_effects, alpha = 0.3, color = NA) +
  theme_classic() +
  # ylim(0,15) +
  labs(x = "Disturbance severity [%]", 
       y = "", title = "p<0.001") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p3


# PlotCaly
p4 <- ggplot(pred_clay, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000),  
              fill = my_color_main_effects, alpha = 0.3, color = NA) +
  theme_classic() +
  scale_y_continuous(breaks = seq(5, 15, 5)) +
  #ylim(0,20) +
  labs(x = "Clay content [%]", y = "", title = "p=0.003") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p4




p_combined_int_no_points <- ggarrange(p1, p4, p2, p3,
                                      labels = c("[a]","[b]", "[c]","[d]"), 
                                      align = 'hv',
                                      font.label = list(size = 8, face = "plain")) # Specify plain font style)

p_combined_int_no_points

# Save the combined plot
ggsave(
  filename = file.path(public_dir, "figs", "Fig2.png"),
  plot = p_combined_int_no_points,
  width = 5.5,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)


#### all predictors for supplement  ------------------------

# update y lables plotting
p.tmp <- p.tmp + labs(y = "")
p1    <- p1 + labs(y = "")
p2    <- p2 + labs(y = "")
p3    <- p3 + labs(y = "")
p4    <- p4 + labs(y = y_lab_reg )


p_combined_int_no_points_supplem <- ggarrange(p.prcp, p.tmp, p1, p4, p2, p3,
                                              labels = c("[a]","[b]", "[c]","[d]", "[e]","[f]"), 
                                              ncol = 3, nrow = 2,
                                              align = 'hv',
                                              font.label = list(size = 8, face = "plain")) # Specify plain font style)

p_combined_int_no_points_supplem

# Save the combined plot
ggsave(
  filename = file.path(public_dir, "figs", "figS_regeneration_by_drivers_supplement.png"),
  plot = p_combined_int_no_points_supplem,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)





# 6. Delayed vs advanced sites: Wilcox plot-----------------------------------------------

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




# 7. Species climate suitability  -------------------------------------------

df_stem_species_class %>% 
  dplyr::filter(stem_density > 0) %>% 
  distinct(Species) %>% 
  count()  #@ 35

# compare current species composition (from field data) with climate suitability from wessely
# Need to adjust tree species naming a bit: 
#   we have less species then wessely: 
# - eg we have besp, which can be betula pendula, betula pubescense
# - same for populus (3), quercus (6), salix ...

# Read species lookup tables
lookup_species_acronyms      <- read.csv(file.path(public_dir, "data", "lookup_species_acronyms.csv"), sep = ';')
lookup_species_field_wessely <- read.csv(file.path(public_dir, "data", "lookup_species_field_wessely.csv"), sep = ';')

# identify what species are present per plot - use simple look up table, 
# add only latin names
species_present_field <- 
  df_stem_species_class %>% 
  dplyr::filter(VegType != "Mature") %>% 
  ungroup() %>% 
  group_by(plot, Species) %>% 
  summarize(reg_stem_density = sum(stem_density, na.rm = T )) %>% 
  mutate(presence = if_else(reg_stem_density > 0, 1, 0)) %>% 
  # dplyr::select(-reg_stem_density ) %>% 
  left_join(lookup_species_acronyms, by = c("Species" = "acc")) %>% 
  dplyr::select(-latin)  %>% 
  dplyr::rename(acc = Species) 

# read table with Wessely species - from species present/absent of plots locations
# only species present on every perios (occurence = 8) are climatically suitable across whole century
species_climate_suitability_raw <- fread(
  file.path(public_dir, "data", "species_presence_clim_change.csv")
)

# add acronyms and consider presence across several species: eg betula sp.
species_suitability_summary <- 
  species_climate_suitability_raw  %>%
  rename(plot = site) %>% 
  # add naming
  left_join(lookup_species_field_wessely, by = c('species' = 'wessely')) %>%
  # group by species to allow occurence of species that have specified genus: eg betula sp.
  group_by(plot, scenario, acc) %>% 
  # Summarize and set sum_presence to 1 if the sum is greater than 1 - 
  # this account for the fact that wessely can have two betulas, I have only 1
  summarize(sum_presence = pmin(sum(overall_presence), 1), .groups = 'drop')

species_suitability_wide <- species_suitability_summary %>%
  pivot_wider(names_from = scenario, values_from = sum_presence) 

# merge both tables: the presently recorded species and species under climate scenarios
df_suitability_observed_vs_modeled <- species_suitability_wide %>% 
  left_join(species_present_field) %>%   # use left join to explude species that are recorded in field, but not present in Wessely database
  dplyr::rename(current = presence)


## Combine field occurence with vs Wessely species databases  ---------
length(unique(df_suitability_observed_vs_modeled$acc))  # final length is 30 species: corss between wessely and observed field database

# create look up table to identify the countries based on regions
region_country_lookup <- df_fin %>%
  group_by(country_pooled) %>%
  dplyr::summarise(region = list(unique(region))) %>%
  unnest(region) %>%
  mutate(region = as.integer(as.character(region)))

# add country indication to merged field and modeled suitability
df_suitability_observed_vs_modeled_reg <- 
  df_suitability_observed_vs_modeled %>%
  # Extract the first two characters of 'plot' as 'region' and convert to integer
  mutate(region = as.integer(substr(plot, 1, 2))) %>%
  # Left join with unique_regions_per_country to get country indication
  left_join(region_country_lookup)



# classify plot level species suitability 
df_suitability_observed_vs_modeled_reg <- df_suitability_observed_vs_modeled_reg %>%
  mutate(
    suitability_rcp26 = case_when(
      current == 1 & rcp26 == 1 ~ "suitable",
      current == 1 & rcp26 == 0 ~ "not_suitable",
      current == 0 & rcp26 == 1 ~ "novel",
      current == 0 & rcp26 == 0 ~ "0"
    ),
    suitability_rcp45 = case_when(
      current == 1 & rcp45 == 1 ~ "suitable",
      current == 1 & rcp45 == 0 ~ "not_suitable",
      current == 0 & rcp45 == 1 ~ "novel",
      current == 0 & rcp45 == 0 ~ "0"
    ),
    suitability_rcp85 = case_when(
      current == 1 & rcp85 == 1 ~ "suitable",
      current == 1 & rcp85 == 0 ~ "not_suitable",
      current == 0 & rcp85 == 1 ~ "novel",
      current == 0 & rcp85 == 0 ~ "0"
    )
  )

###  Overall analysis --------------------------
#### Share of stems that is not suitable under CC scenarios? --------------------------
# Convert suitability columns to logical: TRUE if unsuitable ("0"), FALSE otherwise

anyNA(df_suitability_observed_vs_modeled)  

# get total sum of stems: acdjust for fact that ionly 30 species overlaps between 2 dababases
total_stems         <- sum(df_stem_species_class$stem_density)
total_stems_wessely <- sum(df_suitability_observed_vs_modeled$reg_stem_density )
total_stems_wessely
total_stems
#6172500 - all stems, 5908250 over crossed databases 

# Convert to long format 
df_climate_suitability_regen_species_long <- df_suitability_observed_vs_modeled_reg %>%
  dplyr::select(plot, 
                acc,
                reg_stem_density , 
                country_pooled, 
                suitability_rcp26, 
                suitability_rcp45, 
                suitability_rcp85) %>%
  pivot_longer(cols = starts_with("suitability_rcp"), 
               names_to = "scenario", 
               values_to = "suitability") %>%
  mutate(scenario = gsub("suitability_", "", scenario)) #


df_climate_suitability_regen_species_long %>%  # Remove "suitability_" prefix
  group_by(suitability, scenario) %>% 
  dplyr::filter(suitability == "not_suitable") %>%  
  summarize(sum = sum(reg_stem_density)) %>% 
  mutate(share = round(sum/total_stems*100,1))

# suitability  scenario     sum share
# <chr>        <chr>      <int> <dbl>
#   1 not_suitable rcp26    3782000  61.3
# 2 not_suitable rcp45    4438000  71.9
# 3 not_suitable rcp85    5092750  82.5


#### Most affected species  ---------------
# # Identify species occurence per splots, and which species wuill lose teh most of climate suitability

species_plot_suitability_summary <- 
  df_climate_suitability_regen_species_long %>%
  dplyr::filter(reg_stem_density > 0) %>%
  
  # Count total number of unique plots per species
  group_by(acc) %>%
  mutate(n_plots_current_occurrence  = n_distinct(plot)) %>%
 # View()
  
  # Filter to only suitable cases
  dplyr::filter(suitability == "suitable") %>%
  
  # Count number of suitable plots per species and scenario
  distinct(plot, scenario, acc, n_plots_current_occurrence ) %>%
  count(acc, scenario, name = "n_suitable_plots") %>%
  
  # Join back total plot counts
  left_join(
    df_climate_suitability_regen_species_long %>%
      dplyr::filter(reg_stem_density > 0) %>%
      group_by(acc) %>%
      summarise(n_plots_current_occurrence = n_distinct(plot), .groups = "drop"),
    by = "acc"
  ) %>%
  
  # Calculate share
  mutate(n_unsuitable_plots = n_plots_current_occurrence - n_suitable_plots,
         share_suitable     = round(n_suitable_plots / n_plots_current_occurrence * 100,1),
         share_unsuitable   = 100 - share_suitable, #round(n_unsuitable_plots / n_plots_current_occurrence * 100,1),
         share_overall      = round(n_plots_current_occurrence / n_total_plots * 100,1)) %>%
  #View()
  dplyr::select(scenario, 
                share_overall,
                n_suitable_plots, 
                n_unsuitable_plots,     
                n_plots_current_occurrence,   # number of plots that species occured on
                share_suitable,         
                share_unsuitable,
                acc)

species_plot_suitability_summary

# identify the most affected species (highest share of unsuitable on plot in RCP45)
species_plot_suitability_summary %>%
  dplyr::filter(scenario == "rcp45") %>%
  #dplyr::filter(acc %in% top_species_plot_share) %>% 
  mutate(share_unsuitable = n_unsuitable_plots / n_plots_with_species * 100) %>%
  arrange(desc(share_overall)) %>% # , share_unsuitable
  View()




### Sankey plot test: need to restructure data -----------------

# update naming and coloring
species_colors_grey <- c(
  species_colors,
  other = "#D3D3D3"  # light grey
)

species_labels_grey <- c(
  species_labels,
  other = "Other species"
)


df_top_species_ovr <- 
  df_suitability_observed_vs_modeled_reg %>% 
  dplyr::mutate(
    acc = ifelse(acc %in% top_species_plot_share_v, acc, "other")  # relabel all others
  ) %>%
  #dplyr::filter(acc %in% top_species_plot_share) %>% 
  dplyr::filter(current == 1) %>% 
  dplyr::select(
    plot,
    acc,
    current,
    rcp26,
    rcp45, 
    rcp85,
    reg_stem_density,
    country_pooled
  ) %>% 
  group_by(acc) %>% 
  summarize( n_plots_present = n_distinct(plot),
             avg_stem_dens_current = sum(reg_stem_density)/n_total_plots,
             avg_stem_dens_rcp26 = sum(reg_stem_density[rcp26 == 1]) / n_total_plots,
             avg_stem_dens_rcp45 = sum(reg_stem_density[rcp45 == 1]) / n_total_plots,
             avg_stem_dens_rcp85 = sum(reg_stem_density[rcp85 == 1]) / n_total_plots,
             .groups = "drop") #%>% 


# Prepare long-format data
df_alluvial_stem_density_overall <- 
  df_top_species_ovr %>%
  dplyr::select(acc, 
                avg_stem_dens_current, 
                avg_stem_dens_rcp26, 
                avg_stem_dens_rcp45,
                avg_stem_dens_rcp85) %>%
  pivot_longer(
    cols = starts_with("avg_stem_dens"),
    names_to = "scenario",
    values_to = "stem_density"
  ) %>%
    mutate(
      scenario = case_when(
        scenario == "avg_stem_dens_current" ~ "Current",
        scenario == "avg_stem_dens_rcp26"   ~ "RCP2.6",
        scenario == "avg_stem_dens_rcp45"   ~ "RCP4.5",
        scenario == "avg_stem_dens_rcp85"   ~ "RCP8.5",
        TRUE ~ scenario
      )
    ) %>% 
  mutate(acc = factor(acc, levels = c(top_species_plot_share_v, 'other'))  # Set desired order
  )



# Create the alluvial plot
p_species_stem_density_ovr <- 
  ggplot(df_alluvial_stem_density_overall,
         aes(x = scenario,
             y = stem_density/1000,
             stratum = acc,
             alluvium = acc,
             fill = acc)) +
  geom_flow(alpha = 0.85) +
  geom_stratum() +
  theme_classic2(base_size = 8) +
  labs(title = "",
       x = "Scenario",
       #y =   y_lab_reg, #expression("Reg. stem density [n  1000 " " "*ha^{-1}*"]"),
       y = expression("Regeneration stem density [1000 n "*ha^{-1}*"]"),
       fill = "Species") +
  scale_fill_manual(values = species_colors_grey,
                    labels = species_labels_grey) +
  theme(
    legend.position = "right",
    legend.text = element_text(face = "italic", size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )

p_species_stem_density_ovr


ggsave(
  filename = file.path(public_dir, "figs", "fig_sankey_stem_density_species.png"),
  plot = p_species_stem_density_ovr,
  width = 5,
  height = 3,
  dpi = 300,
  units = "in",
  bg = "white"  # Add white background to ensure clean rendering
)




### Get Tables:  -----------------------------------------------------------
# 1. share of stems/country
# 2. share of plots whene no stemps will remain
# merge into single table

# Summarise stem counts per plot, scenario, country
df_stem_suitability <- 
  df_climate_suitability_regen_species_long %>%
  group_by(country_pooled, scenario) %>%
  summarise(n_plot_country = n_distinct(plot),
    total_stems = sum(reg_stem_density, na.rm = TRUE),
    suitable_stems = sum(reg_stem_density[suitability == "suitable"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    share_suitable = ifelse(total_stems > 0, suitable_stems / total_stems, 0)  # <- Treat 0-stem plots as 0%
  )

# Step 2: Average share across plots (including zero-stem ones as 0%)
df_stem_suitability_avg <- df_stem_suitability %>%
  group_by(country_pooled, scenario, n_plot_country) %>%
  summarise(
    avg_share_suitable = round(mean(share_suitable) * 100, 1),  # As percentage
    .groups = "drop"
  ) %>% 
  mutate(
    country_full = case_when(
      country_pooled == "AT" ~ "Austria",
      country_pooled == "CH" ~ "Switzerland",
      country_pooled == "CZ" ~ "Czech Republic",
      country_pooled == "DE" ~ "Germany",
      country_pooled == "FR" ~ "France",
      country_pooled == "PL" ~ "Poland",
      country_pooled == "SI" ~ "Slovenia",
      country_pooled == "SK" ~ "Slovakia",
      TRUE ~ NA_character_  # Catch anything unexpected
    )
  ) %>%
  arrange(country_full)


df_stem_suitability_avg_wide <- df_stem_suitability_avg %>%
  select(country_pooled, country_full, scenario, avg_share_suitable) %>%
  pivot_wider(
    names_from = scenario,
    values_from = avg_share_suitable,
    names_prefix = "stem_share_"
  ) %>%
  arrange(country_full)


# Result
df_stem_suitability_avg_wide

# TEST END !!!


#### Summary table per country -----------------------------

# get tables: 
# species richness per country 
#### Prepare Table 1: share of climatically suitable tree species per plot/country ------------------

# create a df that contains all plots (also empty ones to calculate properly averages)
df_master <- df_fin %>% 
  dplyr::select(plot, 
                country_pooled )


# compare teh shares of species that are climatically suitable per plot
df_suitability_richness <-   df_suitability_observed_vs_modeled_reg %>% 
  
  group_by(plot, country_pooled) %>%
  mutate(
    
    # calculate number of species
    richness_current = length(unique(acc[current == 1])),
    richness26       = length(unique(acc[current == 1 & rcp26 == 1])),
    richness45       = length(unique(acc[current == 1 & rcp45 == 1])),
    richness85       = length(unique(acc[current == 1 & rcp85 == 1])),
    
    # get shares                     
    rcp26_share = length(unique(acc[current == 1 & rcp26 == 1])) /
      length(unique(acc[current == 1])) * 100,
    rcp45_share = length(unique(acc[current == 1 & rcp45 == 1])) /
      length(unique(acc[current == 1])) * 100,
    rcp85_share = length(unique(acc[current == 1 & rcp85 == 1])) /
      length(unique(acc[current == 1])) * 100
  ) %>%
  ungroup() %>% 
  dplyr::select(plot,   
                acc,
                richness_current,
                richness26,
                richness45,
                richness85,
                current, 
                rcp26, rcp26_share, 
                rcp45, rcp45_share, 
                rcp85, rcp85_share ) %>%
  # filter only present species
  dplyr::filter(current == 1) %>%
  #  View()
  # merge back full table with empty plots
  full_join(df_master) %>% 
  replace_na(list(
    current = 0,
    rcp26_share = 0,
    rcp45_share = 0,
    rcp85_share = 0,
    richness_current = 0,
    richness26 = 0,
    richness45 = 0,
    richness85 = 0
  ))# %>%


#  Calculate shares 
df_suitability_richness_share <- df_suitability_richness %>% #  View()
  mutate(country_pooled = ifelse(is.na(country_pooled), "FR", 
                                 as.character(country_pooled))) %>%
  mutate(country_pooled = as.factor(country_pooled))  %>% # Optional: re-convert to factor if needed
  group_by(country_pooled) %>% 
  #  head()
  summarize(richness   = mean(richness_current),
            richness26 = mean(richness26),
            richness45 = mean(richness45),
            richness85 = mean(richness85),
            richness26_share = mean(rcp26_share),
            richness45_share = mean(rcp45_share),
            richness85_share = mean(rcp85_share)) %>% 
  mutate(
    country_full = case_when(
      country_pooled == "AT" ~ "Austria",
      country_pooled == "CH" ~ "Switzerland",
      country_pooled == "CZ" ~ "Czech Republic",
      country_pooled == "DE" ~ "Germany",
      country_pooled == "FR" ~ "France",
      country_pooled == "PL" ~ "Poland",
      country_pooled == "SI" ~ "Slovenia",
      country_pooled == "SK" ~ "Slovakia",
      TRUE ~ NA_character_  # Catch anything unexpected
    )
  ) %>%
  arrange(country_full)



#### Lost plots:  how many plots per country does not contain any currently present species??? ---------------

df_suitability_plots <- df_suitability_observed_vs_modeled_reg %>%
  group_by(country_pooled, plot) %>%
  mutate(
    any_present = any(current == 1),
    lost26 = if (any_present[1]) all(rcp26[current == 1] == 0) else FALSE,
    lost45 = if (any_present[1]) all(rcp45[current == 1] == 0) else FALSE,
    lost85 = if (any_present[1]) all(rcp85[current == 1] == 0) else FALSE
  ) %>%
  ungroup() %>%
  group_by(country_pooled) %>%
  mutate(
    n_plots_total = n_distinct(plot),
    n_lost_26 = n_distinct(plot[lost26]),
    n_lost_45 = n_distinct(plot[lost45]),
    n_lost_85 = n_distinct(plot[lost85]),
    share_lost_26 = n_lost_26 / n_plots_total * 100,
    share_lost_45 = n_lost_45 / n_plots_total * 100,
    share_lost_85 = n_lost_85 / n_plots_total * 100
  ) %>% 
  dplyr::select(country_pooled,
                n_plots_total,
                n_lost_26,
                n_lost_45,
                n_lost_85,
                share_lost_26,
                share_lost_45,
                share_lost_85) %>% 
  distinct() %>% 
  mutate(
    country_full = case_when(
      country_pooled == "AT" ~ "Austria",
      country_pooled == "CH" ~ "Switzerland",
      country_pooled == "CZ" ~ "Czech Republic",
      country_pooled == "DE" ~ "Germany",
      country_pooled == "FR" ~ "France",
      country_pooled == "PL" ~ "Poland",
      country_pooled == "SI" ~ "Slovenia",
      country_pooled == "SK" ~ "Slovakia",
      TRUE ~ NA_character_  # Catch anything unexpected
    )
  ) %>%
  arrange(country_full)



#### merge climate suitability tables: richenss and lost plots  ----------------------------
df_out <- full_join(df_stem_suitability_avg_wide,
                    df_suitability_plots) %>% 
  full_join(df_suitability_richness_share)

##### simpler version for MS
df_out_MS <- df_out %>% 
  dplyr::select(country_full,
                #n_plots_total,
                stem_share_rcp26 ,
                stem_share_rcp45 , 
                stem_share_rcp85 ,
                share_lost_26,
                share_lost_45,
                share_lost_85)

##### Full version for Supplement 
df_out_supplement <- df_out %>% 
  dplyr::select(-country_pooled) %>% 
  dplyr::select(country_full,
                richness,  # average number of species per plot/country
                richness26,
                stem_share_rcp26 ,
                richness45,
                stem_share_rcp45 , 
                richness85,
                stem_share_rcp85 ,
                n_plots_total,
                n_lost_26,
                share_lost_26,
                n_lost_45,
                share_lost_45,
                n_lost_85,
                share_lost_85)




### Get Summary row: bottom: Totals : values for richness and lost plots across all values:  -----------
#### Richness & remaining stems 
total_richness_stems_full <- df_suitability_observed_vs_modeled_reg %>% 
  summarize(
    # calculate number of species
    richness = length(unique(acc[current == 1])),
    richness26       = length(unique(acc[current == 1 & rcp26 == 1])),
    richness45       = length(unique(acc[current == 1 & rcp45 == 1])),
    richness85       = length(unique(acc[current == 1 & rcp85 == 1])),
    
    # Total stems currently present
    stems_current = sum(reg_stem_density[current == 1], na.rm = TRUE),
    
    # Total stems that are both currently present AND climatically suitable
    stems_suitable_26 = sum(reg_stem_density[current == 1 & rcp26 == 1], na.rm = TRUE),
    stems_suitable_45 = sum(reg_stem_density[current == 1 & rcp45 == 1], na.rm = TRUE),
    stems_suitable_85 = sum(reg_stem_density[current == 1 & rcp85 == 1], na.rm = TRUE)
    
  ) %>% 
  mutate(country_full = "Total",
         richness26_share = richness26/richness*100,
         richness45_share = richness45/richness*100,
         richness85_share = richness85/richness*100,
         stem_share_rcp26 = stems_suitable_26/stems_current*100,
         stem_share_rcp45 = stems_suitable_45/stems_current*100,
         stem_share_rcp85 = stems_suitable_85/stems_current*100
         ) %>% 
  dplyr::select(country_full,
                richness,
                richness26,
                stem_share_rcp26,
                richness45,
                stem_share_rcp45,
                richness85,
                stem_share_rcp85
               )

# make simpler version for MS
total_stem_dens_MS <- total_richness_stems_full %>% 
  dplyr::select(country_full,
                stem_share_rcp26 ,
                stem_share_rcp45 ,
                stem_share_rcp85 )


#### Summary row: Total loss -----------------------
total_lost_full <- 
  df_suitability_observed_vs_modeled_reg %>%
  group_by(plot) %>%
  mutate(
    any_present = any(current == 1),
    lost26 = if (any_present[1]) all(rcp26[current == 1] == 0) else FALSE,
    lost45 = if (any_present[1]) all(rcp45[current == 1] == 0) else FALSE,
    lost85 = if (any_present[1]) all(rcp85[current == 1] == 0) else FALSE
  ) %>%
  ungroup() %>%
  #group_by(country_pooled) %>%
  summarize(
    n_plots_total = n_distinct(plot),
    n_lost_26 = n_distinct(plot[lost26]),
    n_lost_45 = n_distinct(plot[lost45]),
    n_lost_85 = n_distinct(plot[lost85]),
    share_lost_26 = n_lost_26 / n_plots_total * 100,
    share_lost_45 = n_lost_45 / n_plots_total * 100,
    share_lost_85 = n_lost_85 / n_plots_total * 100
  ) %>% 
  dplyr::select(
    n_plots_total,
    n_lost_26,
    share_lost_26,
    n_lost_45,
    share_lost_45,
    n_lost_85,
    share_lost_85) %>% 
  distinct() 


total_lost_MS <- total_lost_full %>% 
  dplyr::select(
    # n_plots_total,
    share_lost_26,
    share_lost_45,
    share_lost_85) 

##### bind total rows for richness and lost plots: for MS and full for supplement 
# total row for MS:  (only percentages)
# full total rown for supplemnt 
total_row_MS         <- cbind(total_richness_MS, total_lost_MS)
total_row_supplement <- cbind(total_richness_full, total_lost_full)


#### Add total rowns to MS and supplement tables: 

# bind the total rows to the original tables: for MS, for Supplement
df_out_MS_with_total         <- bind_rows(df_out_MS, 
                                          total_row_MS) 
df_out_supplement_with_total <- bind_rows(df_out_supplement, 
                                          total_row_supplement)

# format supplement table nicely
df_out_supplement_formatted <- df_out_supplement_with_total %>%
  mutate(
    lost_26 = paste0(n_lost_26, " (", round(share_lost_26, 1), ")"),
    lost_45 = paste0(n_lost_45, " (", round(share_lost_45, 1), ")"),
    lost_85 = paste0(n_lost_85, " (", round(share_lost_85, 1), ")"),
    
    richness = round(richness, 1),
    richness26_fmt = paste0(round(richness26, 1), " (", round(richness26_share, 1), ")"),
    richness45_fmt = paste0(round(richness45, 1), " (", round(richness45_share, 1), ")"),
    richness85_fmt = paste0(round(richness85, 1), " (", round(richness85_share, 1), ")")
  ) %>%
  dplyr::select(
    country_full,
    richness,
    richness26_fmt,
    richness45_fmt,
    richness85_fmt,
    n_plots_total,
    lost_26,
    lost_45,
    lost_85
  )

View(df_out_supplement_formatted)

# sjPlot::tab_df(df_out_MS_with_total,
#                show.rownames = FALSE,
#                file = file.path(public_dir, "tables", "table_MS_clim_suitability_by_country.doc"),
#                digits = 1)
# 
# sjPlot::tab_df(df_out_supplement_formatted,
#                show.rownames = FALSE,
#                file = file.path(public_dir, "tables", "table_Supplement_clim_suitability_by_country.doc"),
#                digits = 1)


# which species will remain?
# Calculate the presence proportion for each scenario
species_presence_proportion <- df_suitability_observed_vs_modeled %>%
  ungroup() %>%
  dplyr::filter(current == 1) %>%
  dplyr::group_by(acc) %>%
  dplyr::summarise(
    current_count = sum(current == 1, na.rm = T),
    rcp26_count = sum(rcp26 == 1, na.rm = T),
    rcp45_count = sum(rcp45 == 1, na.rm = T),
    rcp85_count = sum(rcp85 == 1, na.rm = T),
    current_proportion = (current_count / n_plots) * 100,
    rcp26_proportion = (rcp26_count / n_plots) * 100,
    rcp45_proportion = (rcp45_count / n_plots) * 100,
    rcp85_proportion = (rcp85_count / n_plots) * 100
  ) %>%
  arrange(desc(current_count))

# Display the results
View(species_presence_proportion)


# how may plots have ONLY piab??
plots_with_only_piab <- df_suitability_observed_vs_modeled %>%
  dplyr::filter(reg_stem_density > 0) %>%
  group_by(plot) %>%
  summarise(only_piab = all(acc == "piab"), .groups = "drop") %>%
  dplyr::filter(only_piab) %>%
  pull(plot)  

length(plots_with_only_piab)
# 48 from 849, 5.7%


#### identify teh most affected species --------------------------------
# Evaluate by species: what species are present/ country? and how many of them are not clim suitbale?

# Test for single country: 

current_species<- 
  df_suitability_observed_vs_modeled %>%  
  dplyr::filter(country_pooled == "SK") %>% 
  dplyr::filter(current == 1) %>% 
  #dplyr::filter(rcp26 == 1) %>% 
  ungroup() %>% 
  distinct(acc) #%>%
# rename(current = acc)
#pull()

rcp26_species<- df_suitability_observed_vs_modeled %>%  
  dplyr::filter(country_pooled == "SK") %>% 
  dplyr::filter(rcp26 == 1) %>% 
  ungroup() %>% 
  distinct(acc)  #%>%
#pull() 

# Add an indicator column for presence in each dataset
current_species <- current_species %>% mutate(current_presence = 1)
rcp26_species <- rcp26_species %>% mutate(rcp26_presence = 1)

# Full join to keep all species and display their presence status
merged_species <- full_join(current_species, rcp26_species, by = "acc") %>%
  mutate(current_presence = ifelse(is.na(current_presence), FALSE, current_presence),
         rcp26_presence = ifelse(is.na(rcp26_presence), FALSE, rcp26_presence))

# View the result
View(merged_species)



# Try with cmmapring the vectors:

current_species_v <- 
  df_suitability_observed_vs_modeled %>%  
  dplyr::filter(country_pooled == "SK") %>% 
  dplyr::filter(current == 1) %>% 
  ungroup() %>% 
  distinct(acc) %>%
  pull()

rcp26_species_v<- df_suitability_observed_vs_modeled %>%  
  dplyr::filter(country_pooled == "SK") %>% 
  dplyr::filter(rcp26 == 1) %>% 
  ungroup() %>% 
  distinct(acc)  %>%
  pull() 






current_species_v
rcp26_species_v
length(unique(current_species_v, rcp26_species_v))

species_lost <- setdiff(current_species_v, rcp26_species_v)
length(species_lost)/length(current_species_v)

# Create a data frame that lists all species in the current vector
species_comparison <- data.frame(
  species = current_species_v,
  presence_rcp26 = ifelse(current_species_v %in% rcp26_species_v, "Present", "Lost")
)

# Display the data frame
print(species_comparison)




##### run for all species: -----------------------------------------------------



# Define a function to calculate species loss per climate scenario
calculate_species_loss <- function(data, country, rcp_column) {
  # Filter species present under current conditions
  current_species_v <- data %>%
    dplyr::filter(country_pooled == country, current == 1) %>%
    distinct(acc) %>%
    pull() 
  
  # Filter species present under specified climate scenario
  rcp_species_v <- data %>%
    dplyr::filter(country_pooled == country, .data[[rcp_column]] == 1) %>%
    distinct(acc) %>%
    pull() 
  
  # Calculate species lost and loss share
  species_lost <- setdiff(current_species_v, rcp_species_v)
  loss_share <- length(species_lost) / length(current_species_v)
  n_species_loss <- length(species_lost)
  n_species_current <- length(current_species_v)
  
  # Return results as a list for multiple values
  return(list(loss_share = loss_share, 
              n_species_loss = n_species_loss, 
              n_species_current = n_species_current))
}

# List of climate scenarios to check
rcp_scenarios <- c("rcp26", "rcp45", "rcp85")

# Run the function for each country and scenario, store results in a data frame
species_loss_summary <- df_suitability_observed_vs_modeled %>%
  distinct(country_pooled) %>%
  pull() %>%
  expand.grid(country = ., rcp = rcp_scenarios) %>%
  mutate(rcp = as.character(rcp)) %>%  # Convert rcp to character
  rowwise() %>%
  mutate(
    results = list(calculate_species_loss(df_compare_future_species, country, rcp))
  ) %>%
  unnest_wider(results) %>%  # Separate the list into individual columns
  ungroup()

# Print the resulting summary
print(species_loss_summary)

# Combine `n_species_loss` and `loss_share` into a single formatted column
species_loss_summary <- species_loss_summary %>%
  mutate(
    loss_summary = paste0(n_species_loss, " (", round(loss_share * 100, 1), "%)")
  ) %>%
  dplyr::select(country, n_species_current, rcp, loss_summary) %>% 
  mutate(country = as.character(country)) %>% 
  arrange(country)

# Reshape the data into wide format
species_loss_summary_wide <- species_loss_summary %>%
  pivot_wider(
    names_from = rcp,
    values_from = loss_summary,
    names_prefix = "species_"
  )

# Print the resulting wide-format summary
print(species_loss_summary_wide)

# marge two tables and export

df_species_by_climate <- cbind(species_loss_summary_wide, plot_summary_formatted_share)

df_species_by_climate <- df_species_by_climate %>% 
  dplyr::select(-country_pooled) #%>% 
#dplyr::rename(country = ) #%>% 
#dplyr::rename()


sjPlot::tab_df(df_species_by_climate,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = FALSE,
               file="outTable/climate_suitable_species_plots.doc",
               digits = 1) 

### Barplot of species suitability -------------
# seems wrong ---
df_species_presence <- df_suitability_observed_vs_modeled_reg %>%
  pivot_longer(cols = c(current, rcp26, rcp45, rcp85), 
               names_to = "scenario", 
               values_to = "presence") %>%
  dplyr::filter(presence == 1) %>%  # Keep only species that are present
  dplyr::select(country_pooled, acc, scenario) %>%
  distinct() %>%  # Remove duplicates
  mutate(presence = acc) %>%  # Rename `acc` column for clarity
  pivot_wider(names_from = scenario, values_from = presence, values_fill = "0") 

# make sankey plot:
# how often which species will transit into each climate change scenarios?






# Calculate frequency of occurrence grouped by species, scenario, and presence
df_freq <- df_suitability_observed_vs_modeled_reg %>%
  dplyr::filter(current == 1) %>%
  pivot_longer(cols = c(current, rcp26, rcp45, rcp85), names_to = 'scenario', values_to = 'presence') %>%
  mutate(scenario = factor(scenario, levels = c("current", "rcp26", "rcp45", "rcp85")),
         presence = factor(presence)) %>%
  group_by(acc, scenario, presence) %>%
  summarise(freq = n(), .groups = "drop")



# calculate species richness per country and scanerio - cross link, get a barplot with categories of occurance
# Count unique species in 'current' per country
species_current <- df_suitability_observed_vs_modeled %>%
  dplyr::filter(current == 1) %>%  # Filter where species is present
  group_by(country_pooled) %>%
  dplyr::reframe(n_species_current = unique(acc))

# Count unique species in each future scenario per country
#species_scenarios <- 
df_suitability_observed_vs_modeled %>%
  dplyr::filter(rcp26 == 1 | rcp45 == 1 | rcp85 == 1) %>%  # Keep only species present in any scenario
  group_by(country_pooled) %>%
  dplyr::summarise(
    n_species_rcp26 = unique(acc[rcp26 == 1]),
    n_species_rcp45 = unique(acc[rcp45 == 1]),
    n_species_rcp85 = unique(acc[rcp85 == 1])
  )

# Merge results
species_counts <- species_current %>%
  full_join(species_scenarios, by = "country_pooled") %>%
  arrange(desc(n_species_current))  # Sort by most species currently

# Display results
print(species_counts)










