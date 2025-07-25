# Post-disturbance forest structure

# collect vegetation data from the field
# calculate tree species occurence per plot,
# species richness
# vertical structure
# species stem density 

gc()

# Libs --------------------------------------------------------------------------


library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(stringr)
library(RColorBrewer)
library(scales)
library(colorspace)
library(here)
library(cowplot)
library(sjPlot)


# Input data -------------------------------------------------------------
public_dir <- here("02_analysis/share/EU_forest_regeneration")  # Adjust if needed

source(file.path(public_dir, "code", "00_paths_functions.R"))



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
combined_plot <- cowplot::plot_grid(
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








