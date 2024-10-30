# Analyse the data based on final table: with several speis

# run cluster analysis for climate& environmnet
# for structural data
# get summary tables

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

# library
library(ggridges)

library(purrr)





# stats
library(MASS) #glm.nb for negative binomial models
library(glmmTMB) #Salamanders dataset and lots of families


library(lme4) #pseudo-R2 and and mixed model functionality
library(MuMIn) #dredge function for comparing models, AIC, marginal R^2 for mixed models
library(sjmisc) #pseudo R^2 - but I think gets loaded as a dependency
library(DHARMa) #model diagnostics
library(effects) #what do my marginal effects look like?
library(performance) #binomial model diagnostics
library(emmeans) #post hoc for categorical predictors

#library(AER)  # for hurdle model
library(mgcv)
library(gratia)
library(ggeffects)


# Cluster analysis

library(GGally) # for pairwise comparions of the variables 1:1, also calculates the correlation and sinificance

library(sf)     # read coordinates
library(scales)  # add % for the plots


library(RColorBrewer)


source('my_functions.R')

# Set a global theme -------------------------------------------------------------------

theme_set(
  theme_classic() + 
    theme(
      legend.position = 'bottom',
      text = element_text(size = 4),         # Set all text to size 8
      axis.text = element_text(size = 8),    # Axis tick labels
      axis.title = element_text(size = 8),   # Axis titles
      strip.text = element_text(size = 8),   # Facet labels
      legend.text = element_text(size = 8),  # Legend text
      plot.title = element_text(size = 8)    # Plot title
    )
)



# Read data -----------------------------------------------------------------------

# get vegetation data
load("outData/veg.Rdata")

# final tables on site level
df_fin <- fread('outData/indicators_for_cluster_analysis.csv')

### get coordinates: to add XY coordinates to final model ---------------------
# Replace "your_file.gpkg" with the path to your GPKG file
xy <- st_read("rawData/extracted_soil_data/extracted_soil_data_completed.gpkg")

# Step 1: Create 'site' column by removing last two characters from the 'ID'
#xy <- xy %>%
#  mutate(site = substr(ID, 4, nchar(ID) - 2))

# Step 2: Extract x and y coordinates into separate columns
xy <- xy %>%
  mutate(x = st_coordinates(geom)[, 1],   # Extract x coordinates
         y = st_coordinates(geom)[, 2])   # Extract y coordinates

# Step 3: Group by 'site' and calculate average x and y
xy_site <- xy %>% 
  as.data.frame() %>% 
  rename(site = cluster) %>% 
  group_by(site) %>%
  summarise(x = mean(x), y = mean(y)) %>%
  mutate(site = as.factor(site))
 
dim(xy_site)

# check what is missing (i have already removed the erronerous sites having 3 subplots per plots)
locations          <- unique(xy_site$site)
veg_data_locations <- unique(df_fin$site)

# ad xy coordinates in the final table
df_fin <- df_fin %>% 
  left_join(xy_site, by = join_by(site)) %>% 
  mutate(country_pooled = case_when( country_abbr == "BE" ~ "FR",  # create pooled data for the eco analysis
                                   country_abbr == "LX" ~ "FR",
                                   TRUE~country_abbr))


# Get unique regions for each country_pooled
unique_regions_per_country <- df_fin %>%
  group_by(country_pooled) %>%
  summarise(unique_regions = list(unique(region))) %>% 
  unnest(unique_regions)

# categorize the clim clusters
df_fin <- df_fin %>% 
  mutate(clim_class = case_when(
    clim_cluster_spei3 == 1 ~ "wet-warm-clay",  # Cluster 1: wet, cold, clay
    clim_cluster_spei3 == 2 ~ "hot-dry-sand",  # Cluster 2: hot, dry, sandy (more sand, less clay, more av.nitro than cluster 3)
    clim_cluster_spei3 == 3 ~ "hot-dry-clay"    # Cluster 3: hot, dry, more clay
  )) %>% 
  mutate( clim_cluster_spei3 = as.factor(clim_cluster_spei3),
          clim_class = as.factor(clim_class)) %>% 
  mutate(stem_regeneration = sum_stems_juvenile + sum_stems_sapling)


# export as xy coordinates climate cluster categorized 
df_fin_clim_clust_xy <- st_as_sf(df_fin, coords = c("x", "y"), crs = crs(xy))  

# Step 2: Check the structure of the sf object to ensure everything is correct
print(df_fin_clim_clust_xy)

# Step 3: Export the data as a GeoPackage (GPKG)
#st_write(df_fin_clim_clust_xy, "outData/xy_clim_cluster.gpkg", layer = "df_fin", driver = "GPKG", append = FALSE)

# get only a ataframe of teh climate clusters and sites for easy merging to detiailed veg data
clim_cluster_indicator <- df_fin %>% 
  dplyr::select(site, clim_class)

#fwrite(df_fin, 'outTable/df_fin.csv')

# Analysis ------------------------------------------------------------
# desription of current regeneration state
# structure 
# composition
# 
# What is the current state of teh forest regeneration?
#  stucture? 
#    - stem density
#    - vertical structure
#  Composition? 
#    - richness
#    - species dominance
# H0: presence of regeneration, dominance of teh avanced regeneration
# how often classes cooccur? presence/absence of vertical clases
# higher management intensity, higher species richness
# higher managemnet, less of teh advanced regeneration - 
#   due to the salvage logging and following 

# predictors: investigate, which SPEI is the best? - keep spei12
#             correlate with stem density

# drivers:
# - for stem denisty
# - compare predictors by advanced, middle, delayes

# effect of disturbance size:
#   - higher patch size, less regeneration
#   - higher severity, low regenerations

# environmental condistins;
#  soil: 
#    - more stem density at better soil conditions (higher depth, nutrient content)



## help tables: seral species: --------------

# identify seral stages and whether the species are coniferous or deciduous
df_seral_species <- data.frame(
  Species = c("abal", "lade", "piab", "pist", "pisy", "psme", "taba", "acca", "acpl", "acps", "aehi", "aial", "algl",
              "alin", "alvi", "besp", "cabe", "casa", "fasy", "frex", "fror", "ilaq", "juni", "jure", "osca", "otsp",
              "posp", "potr", "prav", "quro", "qusp", "rops", "saca", "sasp", "soar", "soau", "soto", "tisp", "ulsp"),
  latinn = c("Abies alba", "Larix decidua", "Picea abies", "Pinus strobus", "Pinus sylvestris", "Pseudotsuga menziesii",
             "Taxus baccata", "Acer campestre", "Acer platanoides", "Acer pseudoplatanus", "Aesculus hippocastanum",
             "Ailanthus altissima", "Alnus glutinosa", "Alnus incarna", "Alnus viridis", "Betula spp.", "Carpinus betulus",
             "Castanea sativa", "Fagus sylvatica", "Fraxinus excelsior", "Fraxinus ornus", "Ilex aquifolium", "Juglans nigra",
             "Juglans regia", "Ostrya carpinifolia", "Other species", "Populus spp.", "Populus tremula", "Prunus avium",
             "Quercus robur/petraea", "Quercus spp.", "Robinia pseudoacacia", "Salix caprea", "Salix spp.", "Sorbus aria",
             "Sorbus aucuparia", "Sorbus torminalis", "Tilia spp.", "Ulmus spp."),
  seral_type = c("Late seral", "Early seral", "Late seral", "Late seral", "Early seral", "Late seral", 
                 "Late seral", "Late seral", "Late seral", "Late seral", "Late seral", "Pioneer", "Pioneer",
                 "Pioneer", "Pioneer", "Pioneer", "Late seral", "Late seral", "Late seral", "Pioneer", "Early seral",
                 "Late seral", "Early seral", "Early seral", "Late seral", "Other", "Pioneer", "Pioneer", 
                 "Early seral", "Late seral", "Late seral", "Pioneer", "Pioneer", "Pioneer", "Early seral", 
                 "Early seral", "Early seral", "Late seral", "Late seral"),
  type = c("Coniferous", "Coniferous", "Coniferous", "Coniferous", "Coniferous", "Coniferous", "Coniferous", "Deciduous",
           "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous",
           "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous",
           "Deciduous", "Other", "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous",
           "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous", "Deciduous")
)

# Display the dataframe
df_seral_species

df_stem_sp_sum <- stem_dens_species_long_cluster %>% 
  group_by(cluster, Species) %>% 
  summarise(sum_stem_density = sum(stem_density, na.rm = t)) %>% 
  dplyr::filter(sum_stem_density>0)

# ad indicators for clim cluster
stem_dens_species_long_cluster <- stem_dens_species_long_cluster %>% 
  left_join(clim_cluster_indicator, by = c('cluster' = 'site')) #%>% 


## Species composition overall: total sum of stems ----------------------------------

# Summarize the total stem density per species for each climate class
species_composition_overall <- stem_dens_species_long_cluster %>%
  group_by(Species) %>%
  summarize(sum_stems = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the total stem density per climate class and the share of each species
species_composition_overall <- species_composition_overall %>%
  mutate(total_stems = sum(sum_stems),  # Total stem density in each climate class
         share = (sum_stems / total_stems) * 100) %>%  # Calculate percentage share
  ungroup()


# Find the top 5 species per climate class based on share
top_species_overall <- species_composition_overall %>%
  arrange(desc(share)) %>%  
  slice_head(n = 7) #%>%  # Select the top X species
#dplyr::filter(share > 5)# %>%  # select species with share > 5%

top_species_overall_vect <- top_species_overall$Species

(top_species_overall_vect)


# Species composition by layer -----------------------------------------------------

# Summarize the total stem density per species 
species_composition_layer <- stem_dens_species_long_cluster %>%
  group_by(Species, VegType) %>%
  summarize(sum_stems = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the total stem density per climate class and the share of each species
species_composition_layer <- species_composition_layer %>%
  group_by(VegType) %>%
  mutate(total_stems = sum(sum_stems),  # Total stem density in each climate class
         share = (sum_stems / total_stems) * 100) %>%  # Calculate percentage share
  ungroup()


# Find the top 5 species per climate class based on share
top_species_layer <- species_composition_layer %>%
  group_by(VegType) %>% 
  arrange(desc(share)) %>%  
#  slice_head(n = 7) #%>%  # Select the top X species
dplyr::filter(share > 5)# %>%  # select species with share > 5%

top_species_layer_vect <- top_species_layer$Species

(top_species_overall_vect)
unique((top_species_layer_vect))

# compare teh species by overall ominance and by layer:
setdiff(top_species_layer_vect, top_species_overall_vect )


# Make a barplot - use previous colors and color schemes!!

# Generate custom color palette based on the number of unique species
n_colors <- length(unique(top_species_layer$Species))
my_colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))(n_colors)  # Extend to n_colors

# Manually assign colors to each species based on the desired values
species_colors <- c(
  "frex" = "#A50026",
  "piab" = "#DA362A",
  "pisy" = "#F46D43",
  "fasy" = "#FDAE61",
  "besp" = "#FEE08B",
  "acps" = "#D9EF8B",
  "soau" = "#A6D96A",
  "quro" = "#66BD63",
  "potr" = "#1A9850",
  "abal" = "#006837"
)


# Order the species based on their share
top_species_layer$Species <- reorder(top_species_layer$Species, 
                                     top_species_layer$share, 
                                     decreasing = TRUE)

# Create the bar plot with custom colors and ordered species
p_species_vert_layer <- 
  ggplot(top_species_layer, aes(x = Species, y = share, fill = Species)) +
  geom_bar(stat = "identity", color = 'black') +  # Create bar plot with species share
  facet_grid(VegType ~ ., scales = "free_x") +  # Facet by VegType
  labs(title = "", 
       x = "", y = "Share (%)") +
  scale_fill_manual(values = species_colors) +  # Apply custom color palette
  theme_classic() +  # Use a minimal theme for a clean look
  #coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate the y-axis text (Species)
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Black border around facets
        strip.background = element_rect(color = "black", linewidth = 1),  # Black border around facet labels
        panel.grid.major = element_line(color = "grey", linetype = "dotted"),  # Add grey dashed major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        legend.position = 'none') #+  # Hide the legend
 # Rotate x-axis labels for better readability

### Speceies composition: share per plot: filetr the most important species -------------
# Summarize the total stem density per species 
species_composition_layer_cluster <- stem_dens_species_long_cluster %>%
  group_by(cluster, Species,VegType) %>%
  summarize(sum_stems = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the total stem density per climate class and the share of each species
species_composition_layer_cluster <- species_composition_layer_cluster %>%
  group_by(cluster,VegType) %>%
  mutate(total_stems = sum(sum_stems),  # Total stem density in each climate class
         share = (sum_stems / total_stems) * 100) %>%  # Calculate percentage share
  ungroup() %>% 
  filter(!is.nan(share) & share > 0) 


unique((top_species_layer_vect))

# Reorder Species factor in descending order based on median share
species_composition_layer_cluster <- species_composition_layer_cluster %>%
  dplyr::mutate(Species = factor(Species, levels = names(sort(tapply(share, Species, median), decreasing = TRUE))))



species_composition_layer_cluster %>% 
  dplyr::filter(Species %in% top_species_layer_vect ) %>% 
  ggplot(aes(x = Species, y = share, fill = Species)) +
   stat_summary(fun = median, geom = "bar", color = "black", width = 0.7) +  # Bar plot with median share
 stat_summary(fun.min = function(y) quantile(y, 0.25),  # Lower bound of IQR
              fun.max = function(y) quantile(y, 0.75),  # Upper bound of IQR
              geom = "errorbar", width = 0.2, color = "black") +  # Add IQR error bars
  facet_wrap(~ VegType, scales = "free_x") +  # Facet by VegType
  labs(title = "Median Species Share by Vegetation Type", 
       x = "Species", y = "Median Share (%)") +
  scale_fill_manual(values = species_colors) +  # Apply custom color palette
  theme_classic2() +  # Use a minimal theme for a clean look
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
        legend.position = "none")  # Hide legend if not needed



# get overall number for the species compositions: ----------------------------------

# Calculate median for each species and reorder the factor levels
df_stem_sp_sum_ordered <- df_stem_sp_sum %>%
  dplyr::filter(Species %in% top_species_overall$Species ) %>% 
  dplyr::group_by(Species) %>%
  dplyr::mutate(median_stem_density = median(sum_stem_density, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(Species = reorder(Species, median_stem_density))  # Reorder species by median stem density


# Add a log-transformed column for sum_stem_density
df_stem_sp_sum_ordered <- df_stem_sp_sum_ordered %>%
  mutate(log_sum_stem_density = log10(sum_stem_density + 1))  # Adding 1 to avoid log(0)

# test with log values
df_stem_sp_sum_ordered %>%
  ggplot(aes(x = log_sum_stem_density, y = Species, group = Species)) +
  geom_density_ridges(aes(fill = Species), alpha = 0.5, color = 'NA') +
  scale_fill_manual(values = species_colors) +
  stat_summary(
    aes(x = log_sum_stem_density, fill = Species),  # Add fill aesthetic for inner color
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
    fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
    geom = "pointrange", 
    color = "black",  # Black outline for points
    #fill = "grey",  # Default fill color (or leave as `Species` to map by color)
    shape = 21,  # Shape 21 is a circle with a fill and border
    size = 0.5,
    position = position_nudge(y = 0.5)  # Adjust position slightly
  ) +
  theme_classic() +
  labs(title = "",
       x = "Stem Density (log10)",
       y = "")  +
  theme(legend.position = 'none') +
  scale_x_continuous(
    labels = math_format(10^.x)  # Format x-axis labels as 10^3, 10^4, etc.
  )



# Plot the reordered ridge plot: original values - trimming not working
p_stem_density_species <- df_stem_sp_sum_ordered %>%
  ggplot(aes(x = sum_stem_density, y = Species, group = Species)) +
  geom_density_ridges(aes(fill = Species), alpha = 0.5) +#, quantile_lines = TRUE, quantiles = 2
  stat_summary(
    aes(x = sum_stem_density), 
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
    fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
    geom = "pointrange", 
    color = "black", 
    size = 0.5,
    position = position_nudge(y = .2)  # Adjust position slightly
  ) +
  theme_classic() +
  coord_cartesian(xlim = c(0, 5000)) +
  labs(title = "",
       x = "Stem Density",
       y = "Species")  +
  theme(legend.position = 'none') 

(p_stem_density_species)

ggsave(filename = 'outFigs/p_stem_density_ridge_sum_species.png', 
       plot = p_stem_density_species, 
       width = 4, height = 10, dpi = 300, bg = 'white')



## Plot Stem density per species and vertical class: --------------------------------
# Add a log-transformed column for sum_stem_density
stem_dens_species_long_cluster <- stem_dens_species_long_cluster %>%
  mutate(log_stem_density = log10(stem_density + 1))  # Adding 1 to avoid log(0)

# simplify naming and make new categories for proper stem density
stem_dens_species_long_cluster <- stem_dens_species_long_cluster %>%
  mutate(VegType_acc = case_when(
    VegType == "Mature" ~ "mat",
    VegType == "Juveniles" ~ "juv",
    VegType == "Saplings" ~ "sap",
    TRUE ~ as.character(VegType)  # Keep any other VegType values as they are
  )) %>% 
  mutate(species_VegType = paste(Species, VegType_acc, sep = '_'))


stem_dens_species_long_cluster %>%  
  dplyr::filter(Species %in% top_species_overall_vect[1:5]  ) %>% 
  dplyr::filter(stem_density > 0) %>% 
  mutate(Species = factor(Species, levels = top_species_overall_vect[1:5])) %>%
  ggplot(aes(x = log_stem_density, y = VegType, group = VegType)) +
  geom_density_ridges(aes(fill = VegType), alpha = 0.5) +
  stat_summary(
    aes(x = log_stem_density), 
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
    fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
    geom = "pointrange", 
    color = VegType, 
    size = 0.3,
    position = position_nudge(y = .2)  # Adjust position slightly
  ) +
  facet_grid(Species ~.) +
  theme_classic() +
  coord_cartesian(xlim = c(2,5)) +  # Zoom in on the stem_density range
  labs(
    #title = paste("Density Ridges for Species:", species_name),
    x = "",
    y = ""
  ) +
  theme(legend.position = 'none')


# test density by group???  mutate(Species = factor(Species, levels = top_species_overall_vect[1:5])) %>%
# does not work very well
stem_dens_species_long_cluster %>%  
  dplyr::filter(Species %in% top_species_overall_vect[1:5]  ) %>% 
  dplyr::filter(stem_density > 0) %>% 
  mutate(Species = factor(Species, levels = top_species_overall_vect[1:5])) %>%
  ggplot(aes(x = log_stem_density, group = VegType, fill = VegType)) +#, , y = VegType
  geom_density(alpha=0.6) +
  facet_grid(Species ~.)


# test plotting function ridge density 
plot_ridge_density <- function(data, species_name, xlim_range = c(0, 5)) {
  data %>%
    dplyr::filter(log_stem_density > 0) %>% 
    dplyr::filter(Species == species_name) %>% 
    ggplot(aes(x = log_stem_density, y = VegType, group = VegType)) +
    geom_density_ridges(aes(fill = VegType), alpha = 0.5, trim = T) +
    stat_summary(
      aes(x = log_stem_density), 
      fun = median, 
      fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
      fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
      geom = "pointrange", 
      color = "black", 
      size = 0.5,
      position = position_nudge(y = .2)  # Adjust position slightly
    ) +
    theme_classic() +
    coord_cartesian(xlim = xlim_range) +  # Zoom in on the stem_density range
    labs(
      #title = paste("Density Ridges for Species:", species_name),
      x = "",
      y = ""
    )
}

# Example usage for species 'fasy'
p_piab <- plot_ridge_density(stem_dens_species_long_cluster, species_name = 'piab')
p_fasy <- plot_ridge_density(stem_dens_species_long_cluster, species_name = 'fasy')
p_pisy <- plot_ridge_density(stem_dens_species_long_cluster, species_name = 'pisy')
p_acps <- plot_ridge_density(stem_dens_species_long_cluster, species_name = 'acps')
p_soau <- plot_ridge_density(stem_dens_species_long_cluster, species_name = 'soau')
p_quro <- plot_ridge_density(stem_dens_species_long_cluster, species_name = 'quro')
p_potr <- plot_ridge_density(stem_dens_species_long_cluster, species_name = 'potr')

p_stem_density_ridge <- ggarrange(p_piab,p_fasy,p_pisy,
          p_acps,p_soau,p_quro,p_potr, 
          common.legend = T, ncol = 1,
          labels = top_species_global_vect)#, nrow = 7

p_stem_density_ridge

ggsave(filename = 'outFigs/p_stem_density_ridge.png', 
       plot = p_stem_density_ridge, 
       width = 4, height = 10, dpi = 300, bg = 'white')


## stem density: show only median and IQR?  ----------------------------------


  # Calculate median and IQR for sum_stems per Species and VegType
  df_median_iqr <- stem_dens_species_long_cluster %>%
    dplyr::filter(stem_density >0) %>% 
    group_by(Species, VegType) %>%
    summarize(
      median_stems = median(stem_density, na.rm = TRUE),
      Q1 = quantile(stem_density, 0.25, na.rm = TRUE),  # First quartile (25th percentile)
      Q3 = quantile(stem_density, 0.75, na.rm = TRUE)   # Third quartile (75th percentile)
    ) %>%
    ungroup()
  
  
  
  # Plot the median with IQR for each Species and VegType
  df_median_iqr %>% 
    dplyr::filter(Species %in% c('piab', 'fasy', 'pisy')) %>% 
    ggplot(aes(x = Species, y = median_stems)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.6) +  # Bar plot for median values
    geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, position = position_dodge(0.9)) +  # IQR error bars
    facet_wrap(.~VegType, scales = 'free') +
    labs(title = "Median and IQR of Stem Density per Species and Vegetation Type",
         x = "Species",
         y = "Stem Density (Median ± IQR)") +
    theme_minimal() +
    theme(legend.position = "top")
  
  
## Species composition: vertical class  ------------------------------------
# Summarize the data by VegType and Species
mean_stems_vertical <- stem_dens_species_long_cluster %>%
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

#### Fix the color mapping ----------------------

# Create a color palette based on the number of unique species
n_colors <- length(unique(mean_stems_vertical$Species))
species_colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))(n_colors)

# Map colors to each species
mean_stems_vertical$color <- species_colors[as.numeric(mean_stems_vertical$Species)]



#### Global species composition ----------------------
# Summarize the total stem density per species
species_composition <- stem_dens_species_long_cluster %>%
  group_by(Species) %>%
  summarize(sum_stems = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the total stem density per climate class and the share of each species
species_composition <- species_composition %>%
#  group_by(clim_class) %>%
  mutate(total_stems_clim_class = sum(sum_stems),  # Total stem density in each climate class
         share = (sum_stems / total_stems_clim_class) * 100) %>%  # Calculate percentage share
  ungroup()


#sum of species per fiel wrork: from 37 tree species
# Calculate species richness per clim_class
species_richness <- top_species_per_clim_class %>%
  group_by(clim_class) %>%
  summarise(species_richness = n_distinct(Species))

stem_dens_species_long_cluster %>% 
  ungroup() %>% 
  dplyr::filter(stem_density > 0) %>% 
  summarise(species_richness = n_distinct(Species))
 


# Use different gradient depepnding fof teh seral stage: 

# Find the top 5 species per climate class based on share
top_species_per_clim_class <- species_composition %>%
  #group_by(clim_class) %>%
  arrange(desc(share)) %>%  # Sort species by their share within each climate class
  dplyr::filter(share > 2) %>%  # select species with share > 5%
  left_join(df_seral_species, by = join_by(Species))  #%>% 

# Ensure species are arranged by seral type and within each climate class
top_species_per_clim_class <- top_species_per_clim_class %>%
  arrange(clim_class, seral_type, Species)  # First arrange by seral type and then by Species alphabetically

# Reorder the Species factor based on the seral type
top_species_per_clim_class$Species <- factor(top_species_per_clim_class$Species, 
                                             levels = unique(top_species_per_clim_class$Species[order(top_species_per_clim_class$seral_type)]))



# Load the RColorBrewer package for color palettes
library(RColorBrewer)
# species_table <- table(top_species_per_clim_class$seral_type, top_species_per_clim_class$Species)
# 
# # Count non-zero species in each seral type
# species_counts <- apply(species_table, 1, function(x) sum(x > 0))
# 
# # Extract the count for "Early seral"
# pioneer_n     <- species_counts["Pioneer"]
# early_seral_n <- species_counts["Early seral"]
# late_seral_n  <- species_counts["Late seral"]
# 
# 
# # Define color palettes for each seral type
# pioneer_colors       <- brewer.pal(pioneer_n, "Blues")    # Red shades for pioneer species
# early_seral_colors   <- brewer.pal(early_seral_n, "Oranges")  # Orange shades for early seral species
# late_seral_colors    <- brewer.pal(late_seral_n, "Greens")  # Green shades for late seral species
# 
# # Create a color mapping for Species based on the seral_type
# top_species_per_clim_class$color <- NA
# top_species_per_clim_class$color[top_species_per_clim_class$seral_type == "Pioneer"]     <- pioneer_colors
# top_species_per_clim_class$color[top_species_per_clim_class$seral_type == "Early seral"] <- early_seral_colors
# top_species_per_clim_class$color[top_species_per_clim_class$seral_type == "Late seral"]  <- late_seral_colors
# 
# # Create a named vector for the colors, so each Species has a color
# species_colors <- setNames(top_species_per_clim_class$color, top_species_per_clim_class$Species)

n_colors <- length(unique(top_species_per_clim_class$Species))
my_colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))(n_colors)  # Extend to 12 colors 


# Create the stacked bar plot
p_species_distribution <- ggplot(top_species_per_clim_class, 
                                 aes(x = clim_class, 
                                     y = share, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  geom_text(aes(label = ifelse(share >= 2, paste0(round(share, 1), "%"), "")),
            position = position_stack(vjust = 0.5),  # Labels inside the bars
            size = 3, color = "black") +  # Adjust text size and color
  labs(x = "", y = "Percentage", 
       fill = "Species",
       title = "") +
  scale_fill_manual(values = my_colors) +  # Apply the color palette based on seral type
  theme_classic() +  # Use a clean theme
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

p_species_distribution


ggsave(filename = 'outFigs/fig_p_species_distribution_global.png', 
       plot = p_species_distribution, 
       width = 7, height = 5.5, dpi = 300, bg = 'white')


#### Species compositiosn: country -----------------------------------

# Summarize the total stem density per species for each climate class
species_composition <- stem_dens_species_long_cluster %>%
  group_by(Species, country) %>%
  summarize(sum_stems = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the total stem density per climate class and the share of each species
species_composition <- species_composition %>%
  group_by(country) %>%
  mutate(total_stems_clim_class = sum(sum_stems),  # Total stem density in each climate class
         share = (sum_stems / total_stems_clim_class) * 100) %>%  # Calculate percentage share
  ungroup()


# Find the top 5 species per climate class based on share
top_species_global <- species_composition %>%
 # group_by( country) %>%
  arrange(desc(share)) %>%  # Sort species by their share within each climate class
  slice_head(n = 7) #%>%  # Select the top X species per country
  #dplyr::filter(share > 5)# %>%  # select species with share > 5%

top_species_global_vect <- top_species_global$Species

# Ensure species are arranged by seral type and within each climate class
#top_species_per_clim_class <- top_species_per_clim_class %>%
#  arrange(country, Species)  # First arrange by seral type and then by Species alphabetically

n_colors <- length(unique(top_species_per_clim_class$Species))

my_colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))(n_colors)  # Extend to 12 colors 

# Create the stacked bar plot
p_species_distribution_country <- ggplot(top_species_per_clim_class, 
                                 aes(x = country, 
                                     y = share, 
                                     fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  geom_text(aes(label = ifelse(share >= 2, paste0(round(share, 1), "%"), "")),
            position = position_stack(vjust = 0.5),  # Labels inside the bars
            size = 3, color = "black") +  # Adjust text size and color
  labs(x = "", y = "Percentage", 
       fill = "Species",
       title = "") +
  scale_fill_manual(values = my_colors) +  # Apply the color palette based on seral type
  theme_classic() +  # Use a clean theme
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) 

p_species_distribution_country


ggsave(filename = 'outFigs/fig_p_species_distribution_global_country.png', 
       plot = p_species_distribution_country, 
       width = 7, height = 5, dpi = 300, bg = 'white')




library(treemapify)

stem_dens_species_long_cluster %>% 
  group_by(Species, VegType) %>% 
  mutate(VegType = factor(VegType, levels = c('Saplings', "Juveniles", "Mature"))) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(VegType) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share       = sum_stems/sum_vegType*100) %>%
  dplyr::filter(share>5) %>% 
  #slice_max(order_by = share, n = 10) %>% 
  ggplot(aes(fill = Species, 
             area = sum_stems ,
             label = paste(Species, "\n", round(share,0)))) +
  geom_treemap() +
  geom_treemap_text(colour ="white", place = "centre") +
  facet_grid(~VegType)
# geom_mosaic()
#geom_mosaic(aes(x = product(VegType), 
#                fill = do_you_recline))

## Structure ---------------------------------------------------------

### Vertical structure  --------------------------------------

#### Get presence absence data for indiviual layers --------------------------

# Group by cluster and VegType, check if stem_density > 0
vert_class_presence_absence <- stem_dens_species_long_cluster %>%
  group_by(cluster, VegType) %>%
  summarise(has_stem_density = sum(stem_density > 0, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  # Filter only where there is stem density present
  dplyr::filter(has_stem_density) %>%
  dplyr::select(cluster, VegType) %>% 
  mutate(Presence = 1) %>%  # Add a column with 1 indicating presence
  pivot_wider(names_from = VegType, values_from = Presence, values_fill = 0)  # Convert to wide format, filling NAs with 0



# Now find clusters where no VegType has stem_density > 0
# Find clusters where no vegType has stem_density
all_clusters <- stem_dens_species_long_cluster %>%
  ungroup() %>% 
  dplyr::select(cluster, VegType) %>% 
  distinct(cluster)

# Combine clusters with and without stem density
vert_class_presence_absence_fin <- vert_class_presence_absence  %>% 
  right_join(all_clusters)



# visualize vertical classes by UpSEt plot

# Upset plot -----------------------------------

library(UpSetR)



# full data:

# Step 2: Select only the columns with presence/absence data
upset_data <- vert_class_presence_absence_fin %>% dplyr::select(Saplings, Juveniles, Mature) %>% 
  as.data.frame()

# Step 3: Create the UpSet plot
upset(upset_data, sets = c('Saplings', 'Juveniles', 'Mature'), order.by = "freq")


# upset data test ---------------------------------------------------------

dd <- data.frame(site = c(1,2,2,3,3,3,4,5,5,6,7,7,7),
                 vert = c('m', 
                          's', 'j',
                          'j','s','m',
                          's',
                          'j','s',
                          's',
                          'm','j','s'))
# Step 1: Create a binary presence/absence matrix for each site
dd_wide <- dd %>%
  pivot_wider(names_from = vert, values_from = vert, 
              values_fn = length, values_fill = 0) %>%
  mutate(m = ifelse(m > 0, 1, 0),
         j = ifelse(j > 0, 1, 0),
         s = ifelse(s > 0, 1, 0))

# Step 2: Select only the columns with presence/absence data
upset_data <- dd_wide %>% dplyr::select(m, j, s) %>% 
  as.data.frame()

# Step 3: Create the UpSet plot
upset(upset_data, sets = c("m", "j", "s"), order.by = "freq")






# PLOT: summary inicators -----------------------------------------------------
# Define the function to create a customizable violin plot
create_violin_plot <- function(df, y_var, y_label) {
  df %>%
    ggplot(aes(x = clim_class, y = !!sym(y_var))) +  # Use dynamic y variable
    geom_violin(aes(fill = clim_class), 
                trim = TRUE, alpha = 0.6) +  # Create the violin plot
    geom_boxplot(width = 0.15) +
    scale_fill_manual(values = c("orange", "red", "blue")) +
    labs(x = '',
         fill = 'Clim_ENV cluster',
         y = y_label) +  # Use dynamic y label
    theme_classic2() +
    theme(legend.position = "none")
}

df_fin$plot_stem_density <- df_fin$stem_density/1000
df_fin$plot_rIVI <- df_fin$rIVI*100
# Example usage for rIVI, n_vertical, and richness

# For rIVI
p_viol_stem_density <- create_violin_plot(df_fin, "plot_stem_density", "Stem density [#*1000]")

# For rIVI
p_viol_rIVI <- create_violin_plot(df_fin, "plot_rIVI", "rIVI [%]")

# For n_vertical
p_viol_vert <- df_fin %>%
  ggplot(aes(x = clim_class, y = n_vertical)) +  # Use dynamic y variable
  geom_violin(aes(fill = clim_class), 
              trim = TRUE, alpha = 0.6) +  # Create the violin plot
  # geom_boxplot(width = 0.15) +
  scale_fill_manual(values = c("orange", "red", "blue")) +
  labs(x = '',
       fill = 'Clim_ENV cluster',
       y = 'Vertical layers [#]') + 
  theme_classic2() +
  theme(legend.position = "none") #create_violin_plot(df_fin, "n_vertical", "Vertical layers [#]")

# For richness
p_viol_richness <- create_violin_plot(df_fin, "richness", "Richness [#]")

p_indicators_violin <- ggarrange(p_viol_stem_density,p_viol_vert, p_viol_rIVI, p_viol_richness,
                                 labels = c("[a]", "[b]", "[c]", "[d]"),  # Add labels for each subplot,
                                 legend = "none",
                                 font.label = list(size = 10, face = "plain") ) # Adjust legend position if needed)

(p_indicators_violin)
ggsave(filename = 'outFigs/fig_indicators_violin.png', 
       plot = p_indicators_violin, 
       width = 7, height = 7, dpi = 300, bg = 'white')





## Descriptive plots  -----------------------------------------------------------

# Richness desc 

# decrsibe richness: how many clusters have how many species?
df_fin %>% 
  group_by(richness) %>% 
  dplyr::summarize(count = n(),
                   prop = count/cluster_n*100)


# what are dominant species??
df_fin %>% 
  dplyr::filter(rIVI> 0.5) %>%
  # nrow()
  group_by(dominant_species) %>% 
  dplyr::summarize(count = n(),
                   prop = count/cluster_n*100) %>% 
  arrange(desc(prop)) %>% 
  mutate(dominant_species = factor(dominant_species, 
                                   levels = unique(dominant_species[order(count, decreasing = TRUE)]))) %>%
  ggplot(aes(x = dominant_species,
             y = count)) +
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = paste0(round(prop, 1), "%")), 
            position = position_stack(vjust = 1.4), 
            angle = 90,
            color = "black",
            size= 4) + # Adjust color as needed
  labs(x = '',
       y = 'Counts [#]') +
  theme_bw()




# 3. how many clusters/plots has less then 2000 regen/ha? ----------------------------------
# which country is the most affected? (where they ccur most often?)



df_stock_density <- 
  df_stems %>% 
  ungroup(.) %>% 
  #  dplyr::filter(cluster %in% cluster_id_keep) %>% 
  mutate(sum_stems = as.integer(sum_stems),
         dens_category = case_when(sum_stems == 0 ~ '0',
                                   #sum_stems == 500 ~ '1-500',
                                   #sum_stems == 1000 ~ '500-1000',
                                   #sum_stems == 1500 ~ '1500-2000',
                                   sum_stems >= 1   & sum_stems < 501 ~ '1-500',
                                   sum_stems >= 501 & sum_stems < 1001 ~ '500-1000',
                                   sum_stems >= 1001 & sum_stems < 1501 ~ '1000-1500',
                                   sum_stems >= 1501 & sum_stems < 2001 ~ '1500-2000',
                                   sum_stems >= 2001 ~ '>2000'
         ),
         dens_class = case_when(sum_stems < 2001 ~ 'low',
                                sum_stems >= 2001 ~ 'sufficient'
         )) %>%
  mutate(dens_category = factor(dens_category, levels = 
                                  c('0', '1-500','500-1000','1000-1500', '1500-2000', '>2000' )),
         dens_class = factor(dens_class, levels = c('low', 'sufficient')))# %>% 
#View(df_stock_density)

# make a table for potential regeneration delay - link with climate data
# make a barplots
df_regen_delay <- df_fin %>% 
  mutate(reg_delay_class = case_when(stem_density  < 2001 ~ 'low',
                                     stem_density  >= 2001 ~ 'sufficient')) %>% 
  mutate(reg_delay_class = factor(reg_delay_class, levels = c('low', 'sufficient')))


median_iqr <- function(y) {
  data.frame(
    y = median(y, na.rm = TRUE),
    ymin = quantile(y, probs = 0.25, na.rm = TRUE),
    ymax = quantile(y, probs = 0.75, na.rm = TRUE)
  )
}



# Define the function
create_plot <- function(data, 
                        x_var, 
                        y_var, 
                        #x_label = "", 
                        y_label = y_var,
                        x_annotate = 0, lab_annotate = "lab ann") {
  
  p <- ggplot(data, aes_string(x = x_var, #"reg_delay_class", 
                               y = y_var)) +
    stat_summary(fun.data = median_iqr, geom = "pointrange") +
    stat_summary(fun = median, geom = "point") +
    labs(xlab = NULL) +
    annotate("text", x = x_annotate, y = Inf, label = lab_annotate, hjust = 0.5, vjust = 1.5) +
    theme_minimal(base_size = 10) +
    theme(aspect.ratio = 1, 
          #legend.position = 'bottom',
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          #axis.text.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 8)) 
  
  return(p)
}





## for climate ------------------------------------------------------------------
# Using the function to create plots
p1 <- create_plot(df_regen_delay, x_var = "reg_delay_class", 
                  y_var = "tmp",  x_annotate = 1.5, lab_annotate = "0.06")
p2 <- create_plot(df_regen_delay, 
                  x_var = "reg_delay_class", 
                  y_var = "prec",  x_annotate = 1.5, lab_annotate = "***")
p3 <- create_plot(df_regen_delay, 
                  x_var = "reg_delay_class", 
                  y_var = "spei", x_annotate = 1.5, lab_annotate = "***")

windows(7,3)
ggarrange(p1, p2, p3, ncol = 3, labels = c('[a]', '[b]', '[c]'))

# get summary statistics
df_regen_delay %>% 
  ungroup() %>% 
  group_by(reg_delay_class) %>% 
  summarize(
    spei_med = median(spei, na.rm = TRUE),
    spei_sd = sd(spei, na.rm = TRUE),
    spei_25 = quantile(spei, 0.25, na.rm = TRUE), # rich_mean -rich_sd, #
    spei_75 = quantile(spei, 0.75, na.rm = TRUE), # rich_mean +rich_sd,
    tmp_med = median(tmp   , na.rm = TRUE),
    tmp_sd   = sd(tmp  , na.rm = TRUE),
    tmp_25 = quantile(tmp  , 0.25, na.rm = TRUE), # dens_mean - dens_sd, #
    tmp_75 = quantile(tmp  , 0.75, na.rm = TRUE), # dens_mean + dens_sd, #
    prec_med = median(prec   , na.rm = TRUE),
    prec_sd   = sd(prec  , na.rm = TRUE),
    prec_25 = quantile(prec  , 0.25, na.rm = TRUE), # dens_mean - dens_sd, #
    prec_75 = quantile(prec  , 0.75, na.rm = TRUE), 
    .groups = 'drop')


## test differences between medians of teh classes -----------------------------

##### Loop over wilcx test

perform_wilcox_test <- function(data, variable_name, group1, group2) {
  # Subset variable values for each specified 'reg_delay_class' group
  group1_values <- data[[variable_name]][data$reg_delay_class == group1]
  group2_values <- data[[variable_name]][data$reg_delay_class == group2]
  
  # Perform the Mann-Whitney U test between the two subsets
  test_result <- wilcox.test(group1_values, group2_values)
  
  # Return the test result
  return(test_result)
}

# Example usage for 'spei' comparing 'low' vs 'sufficient'
test_result_spei <- perform_wilcox_test(df_regen_delay, "spei", "low", "sufficient")
print(test_result_spei)

# Example usage for 'temp' comparing 'low' vs 'sufficient'
test_result_temp <- perform_wilcox_test(df_regen_delay, "tmp", "low", "sufficient")
print(test_result_temp)

# Example usage for 'prcp' comparing 'low' vs 'sufficient'
test_result_prcp <- perform_wilcox_test(df_regen_delay, "prec", "low", "sufficient")
print(test_result_prcp)




## For disturbance size and intensity -------------------------------------------

# Using the function to create plots
p1 <- create_plot(df_regen_delay, x_var = "reg_delay_class", 
                  y_var = "disturbance_severity",  x_annotate = 1.5, lab_annotate = "n.s.")
p2 <- create_plot(df_regen_delay, 
                  x_var = "reg_delay_class", 
                  y_var = "distance_edge",  x_annotate = 1.5, lab_annotate = "n.s.", y_label = 'distance edge [m]')

windows(5,3)
ggarrange(p1, p2,  ncol = 2, labels = c('[a]', '[b]'))

hist(df_regen_delay$disturbance_severity)

# get summary statistics
df_regen_delay %>% 
  ungroup() %>% 
  group_by(reg_delay_class) %>% 
  summarize(
    severity_med = median(disturbance_severity, na.rm = TRUE),
    severity_sd = sd(disturbance_severity, na.rm = TRUE),
    severity_25 = quantile(disturbance_severity, 0.25, na.rm = TRUE), # rich_mean -rich_sd, #
    severity_75 = quantile(disturbance_severity, 0.75, na.rm = TRUE), # rich_mean +rich_sd,
    dist_edge_med = median(distance_edge   , na.rm = TRUE),
    dist_edge_sd   = sd(distance_edge  , na.rm = TRUE),
    dist_edge_25 = quantile(distance_edge  , 0.25, na.rm = TRUE), 
    dist_edge_75 = quantile(distance_edge  , 0.75, na.rm = TRUE), 
    .groups = 'drop')


## test differences between medians of teh classes -----------------------------

##### Loop over wilcx test

perform_wilcox_test <- function(data, variable_name, group1, group2) {
  # Subset variable values for each specified 'reg_delay_class' group
  group1_values <- data[[variable_name]][data$reg_delay_class == group1]
  group2_values <- data[[variable_name]][data$reg_delay_class == group2]
  
  # Perform the Mann-Whitney U test between the two subsets
  test_result <- wilcox.test(group1_values, group2_values)
  
  # Return the test result
  return(test_result)
}

# Example usage for 'spei' comparing 'low' vs 'sufficient'
perform_wilcox_test(df_regen_delay, "disturbance_severity", "low", "sufficient")

# Example usage for 'temp' comparing 'low' vs 'sufficient'
perform_wilcox_test(df_regen_delay, "distance_edge", "low", "sufficient")



# make a plot? not working wnow
df_stock_density %>% 
  ggplot(aes(#x = manag,
    # y = dens_category,
    fill = dens_category)) +
  geom_bar(position = 'fill') +
  scale_fill_manual(values = c("0" = "red", 
                               "500-1000" = "lightgreen", 
                               "1000-2000" = "green", 
                               ">2000" = "darkgreen",
                               "1-500" = "forestgreen")) +
  ylab('Share [%]') + 
  theme_classic() +
  labs(fill = "Density class") +
  xlab('')

unique(df_stock_density$dens_category)

df_stock_density %>% 
  group_by(dens_category) %>% 
  summarise(n = n()) %>%
  mutate(share = round(n/n_clusters*100,1)) %>% 
  ggplot(aes(fill = dens_category     , 
             area = share      ,
             label = paste(dens_category, "\n", round(share,2), "%"))) +
  geom_treemap() +
  geom_treemap_text(colour ="white", place = "centre")# +
#facet_grid(~VegType)




#### Pearson chi squared test --------------------------
# Create a contingency table
contingency_table <- table(df_stock_density$manag, df_stock_density$dens_category)

# Perform the chi-squared test
chi_squared_test <- chisq.test(contingency_table)

# View the results of the chi-squared test
print(chi_squared_test)

# invalid chi suqre test - 0 categories, very little counts (less then 5)
# Combine sparse categories
df_stock_density$combined_category <- ifelse(df_stock_density$dens_category %in% c("0", "1-500"), 
                                             "0-1-500", 
                                             df_stock_density$dens_category)

# Create a new contingency table with combined categories
new_contingency_table <- table(df_stock_density$manag, df_stock_density$combined_category)

# Perform the chi-squared test on the new table
new_chi_squared_test <- chisq.test(new_contingency_table)

# View the results
print(new_chi_squared_test)



# Calculate the percentages
df_stock_density_simpl <- df_stock_density %>%
  ungroup(.) %>% 
  count(manag, dens_category) %>%
  group_by(manag) %>%
  mutate(freq = n / sum(n),
         label = scales::percent(freq, accuracy = 0.1))


df_stock_density_simpl %>% 
  ggplot( aes(x = manag, y = freq, fill = dens_category)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = label, y = freq/2), position = position_fill(vjust = 0.5), size = 4) +
  ylab('Share [%]') +
  theme_classic() +
  labs(fill = "Density class")





# Compare categories occurences by chi square






# 4. stems vs weather - mean 2018-2023 ---------------------------------

df_fin
windows()
pairs(stem_density    ~ tmp + tmp_z + prec + prcp_z + management_intensity+ salvage_intensity + protection_intensity + distance_edge, df_fin)

pairs(stem_density    ~  management_intensity+ salvage_intensity + protection_intensity, df_fin)


# test corelations 

plot(df_fin$tmp, df_fin$spei)
plot(df_fin$prec, df_fin$spei)
plot(df_fin$tmp, df_fin$prec)

plot( df_fin$management_intensity, df_fin$stem_density)
plot( df_fin$salvage_intensity, df_fin$stem_density)
plot( df_fin$protection_intensity, df_fin$stem_density)


# keep only spei instead of the tmp and prec

# keep also annomalies?
plot(df_fin$tmp_z, df_fin$spei)

# test, which one of teh variables exaplin teh stem density better?




cor(df_fin$tmp_z, df_fin$spei, method = 'spearman')
cor(df_fin$tmp, df_fin$spei, method = 'spearman')
#[1] -0.7187752
cor(df_fin$prec, df_fin$spei, method = 'spearman')
# [1] 0.2349638

# check how management is correlated?
cor(df_fin$management_intensity, df_fin$salvage_intensity, method = 'spearman')
cor(df_fin$protection_intensity, df_fin$salvage_intensity, method = 'spearman')



# decide: which parameters are better: temp_z, temp, prcp, prcp_z or spei?? or their interation?

# Calculate correlation matrix for predictors
correlation_matrix <- cor(df_fin[, c('tmp_z', 'tmp', 'prec', 'prcp_z', 'spei')])


# remove spei
correlation_matrix <- cor(df_fin[, c('tmp_z', 'tmp', 'prec', 'prcp_z')], method = 'spearman')

correlation_matrix




#5. Drivers ---------------------------------
# for all regeneration (juveniles and saplings)
# split table in two: drivers for advanced (> 1000 stems of juveniles/ha)
#                     drivers for delayed regeneration (<50 stems/ha: saplings  + juveniles)  

## Prep final table --------------

df_fin <- df_fin %>% 
  mutate(country_full    = factor(country_full),
         country_abbr    = factor(country_abbr),
         country_pooled  = factor(country_pooled),
         region          = factor(region))



df_fin <- df_fin %>% 
  mutate(
    tmp_c = tmp - mean(tmp, na.rm = TRUE),  # Centering temperature
    prcp_c = prcp - mean(prcp, na.rm = TRUE) # Centering precipitation
  )

# account for climatic variables - aggregated on 9 km resolution:
df_fin <- df_fin %>%
  mutate(clim_grid = factor(paste(tmp, prcp, sep = "_")))  # Create a unique identifier for temp/precip combinations



# Make sure df_fin has no missing values, if necessary
df_fin <- na.omit(df_fin)

# make  regeneration classes: advanced vs delayed
df_fin <- df_fin %>% 
  mutate(adv_delayed = factor(ifelse(stem_regeneration <= 50, "Delayed", 
                                     ifelse(sum_stems_juvenile >= 1000, "Advanced", "Other")),
                              levels = c("Delayed", "Other", "Advanced"))) %>%  
  # make larger categories to differentiate between groups
  mutate(adv_delayed_wider = factor(ifelse(stem_regeneration <= 500, "Delayed", 
                                     ifelse(sum_stems_juvenile >= 1000 , "Advanced", #| stem_regeneration >= 2000
                                            "Other")),
                              levels = c("Delayed", "Other", "Advanced"))) %>%  
  
  # create binary classes for advanced vs delayed
  mutate(delayed = ifelse(stem_regeneration <= 50, 1, 0),
         advanced = ifelse(sum_stems_juvenile >=  1000, 1, 0))

# check my categories?
df_fin %>% 
  dplyr::select(site, 
                stem_density,
                stem_regeneration,
                sum_stems_juvenile,
                sum_stems_sapling, 
                #sum_stems_mature,
                adv_delayed,
                adv_delayed_wider
  ) %>% 
  View()




## Variables selection: Correlations  ------------------------------------------------------

# Select relevant columns for the plot
df_pairs <- df_fin %>%
  dplyr::select(stem_density, spei1, spei3, spei6, spei12, spei24,
                drought_spei1, drought_spei3, drought_spei6, drought_spei12, drought_spei24)

# Create a pairs plot with correlations, distributions, and scatter plots
ggpairs(df_pairs,
        lower = list(continuous = "smooth"),  # Scatter plots with smoothing lines
        diag = list(continuous = "barDiag"),  # Histograms on the diagonal
        upper = list(continuous = "cor"),     # Correlation coefficients in the upper triangle
        title = "Pairs Plot: Stem Density and SPEI Scales"
)


# very little correlations: try spearman

# Calculate Spearman correlations between stem_density and each SPEI scale
spearman_correlations <- df_fin %>%
  dplyr::select(stem_density, spei1, spei3, spei6, spei12, spei24,
                drought_spei1, drought_spei3, drought_spei6, drought_spei12, drought_spei24,
                tmp, tmp_z, prcp, prcp_z,
                drought_tmp, drought_prcp) %>%
  summarise(
    sp_spei1 = cor(stem_density, spei1, method = "spearman", use = "complete.obs"),
    sp_spei3 = cor(stem_density, spei3, method = "spearman", use = "complete.obs"),
    sp_spei6 = cor(stem_density, spei6, method = "spearman", use = "complete.obs"),
    sp_spei12 = cor(stem_density, spei12, method = "spearman", use = "complete.obs"),
    sp_spei24 = cor(stem_density, spei24, method = "spearman", use = "complete.obs"),
    sp_drought_spei1 = cor(stem_density, drought_spei1, method = "spearman", use = "complete.obs"),
    sp_drought_spei3 = cor(stem_density, drought_spei3, method = "spearman", use = "complete.obs"),
    sp_drought_spei6 = cor(stem_density, drought_spei6, method = "spearman", use = "complete.obs"),
    sp_drought_spei12 = cor(stem_density, drought_spei12, method = "spearman", use = "complete.obs"),
    sp_drought_spei24 = cor(stem_density, drought_spei24, method = "spearman", use = "complete.obs"),
    sp_tmp = cor(stem_density, tmp, method = "spearman", use = "complete.obs"),
    sp_tmp_z = cor(stem_density, tmp_z, method = "spearman", use = "complete.obs"),
    sp_prcp = cor(stem_density, prcp, method = "spearman", use = "complete.obs"),
    sp_prcp_z = cor(stem_density, prcp_z, method = "spearman", use = "complete.obs"),
    sp_drought_tmp = cor(stem_density, drought_tmp, method = "spearman", use = "complete.obs"),
    sp_drought_prcp = cor(stem_density, drought_prcp, method = "spearman", use = "complete.obs"),
  )

# Print the Spearman correlations
# Transpose the data frame so that each variable becomes a row
spearman_transposed <- t(spearman_correlations)

# Convert to a data frame and add column names
spearman_transposed <- data.frame(value = spearman_transposed[, 1], row.names = rownames(spearman_transposed))

# Order the transposed values
ordered_spearman <- spearman_transposed[order(spearman_transposed$value), , drop = FALSE]

# Print the ordered values and their corresponding variables
print(ordered_spearman)

# precipitation is the most correlated to stem_density

# > print(ordered_spearman)
# value
# sp_drought_tmp    -0.035620505
# sp_tmp            -0.021468492
# sp_drought_spei24 -0.016574555
# sp_tmp_z           0.002567999
# sp_drought_spei12  0.049451497
# sp_spei1           0.083722968
# sp_spei24          0.092396768
# sp_spei12          0.096146369
# sp_prcp_z          0.096295058
# sp_spei3           0.105548424
# sp_spei6           0.111553402
# sp_drought_spei1   0.118807397
# sp_drought_spei3   0.124401148
# sp_drought_spei6   0.162161276
# sp_drought_prcp    0.203395486
# sp_prcp            0.210360501


#### check for multicollinearity -----------------------------------------------------

library(car)
selected_data <- df_fin %>%
  dplyr::select(stem_density, # spei1, spei3, spei6, spei12, spei24,
                #drought_spei1, drought_spei3, 
                #                drought_spei6, 
                #drought_spei12, drought_spei24,
                tmp, #tmp_z, 
                prcp #, #prcp_z,
                # drought_tmp, 
                # drought_prcp
  )

# Step 2: Fit a linear model predicting stem_density
lm_model <- lm(stem_density ~ ., data = selected_data)

# Step 3: Run VIF to check multicollinearity
vif_values <- vif(lm_model)

(vif_values)

# final climatic predictors are tmp and prec

### test with different families ----------------------------------------------
hist(df_fin$stem_density)


# Fit a GAM model with a Negative Binomial distribution
m1 <- gam(stem_density ~ s(spei6, k = 15),  # Factors included without s() for categorical variables
          family = nb,  # Negative Binomial to handle overdispersion
          data = df_fin)



# Fit a GAM with Tweedie distribution (useful for zero-inflation)
m.tw.juv <- gam(sum_stems_juvenile ~ s(drought_spei6, k = 10), 
                family = tw,  # Adjust var.power based on data
                data = df_fin)

m.tw.mature <- gam(sum_stems_mature ~ s(drought_spei12, k = 10), 
                   family = tw,  # Adjust var.power based on data
                   data = df_fin)

m.tw.sapl <- gam(sum_stems_sapling ~ s(drought_spei3, k = 10), 
                 family = tw,  # Adjust var.power based on data
                 data = df_fin)


hist(df_fin$sum_stems_mature)

AIC(m1, m.tw1)

appraise(m.tw.mature)
summary(m.tw.mature)
draw(m.tw.mature)
gam.check(m.tw1)
k.check(m.tw1)

# TW has a better fit, also can handle zero!

##### Run univariate models for set of dependent variables (stem density) ------------------

# List of dependent variables
dependent_vars <- c("sum_stems_juvenile", 
                    "sum_stems_sapling", 
                    "stem_regeneration", # sum of juveniles and saplings
                    #"sum_stems_mature",
                    "advanced",
                    "delayed"
                    #"stem_density"
)

# List of predictor variables (spei1 to spei24, drought_spei1 to drought_spei24)
predictor_vars <- c("spei1", "spei3", "spei6", 
  "spei12", 
  #"spei24",
                    "tmp", 
                    #"tmp_z", 
                    "prcp", 
                    #"prcp_z", 
                   # "drought_spei1", "drought_spei3", "drought_spei6", 
  #                  "drought_spei12", 
  "drought_spei24",
   #                 "drought_tmp", "drought_prcp",
                    "management_intensity",
                    "distance_edge", 
                    "disturbance_severity", 
                    "clay_extract", 
                    "clay_extract", 
                    "depth_extract", 
                    "av.nitro")


# Initialize a data frame to store AIC values and deviance explained
model_metrics <- data.frame(Predictor = character(), 
                            Dependent = character(), 
                            AIC = numeric(), DevianceExplained = numeric())



# Loop over each dependent variable
for (dep in dependent_vars) {
  #print(dep)
  # Loop over each predictor
  for (pred in predictor_vars) {
    #print(pred)
    # Fit the model
    formula <- as.formula(paste(dep, "~ s(", pred, ", k = 10)"))
    #print(formula)
    model <- gam(formula, family = tw(), method = 'REML', data = df_fin)
    
    # Extract model summary
    model_summary <- summary(model)
    
    # Store the AIC value and deviance explained
    model_metrics <- rbind(model_metrics, 
                           data.frame(Predictor = pred, 
                                      Dependent = dep, 
                                      AIC = AIC(model), 
                                      DevianceExplained = round(model_summary$dev.expl*100,1)))
  }
}

# View the AIC values and deviance explained
View(model_metrics)

# # Select the best predictor for each dependent variable based on the lowest AIC
best_predictors <- 
  model_metrics %>%
  mutate(dependent_category = case_when(
    grepl("sapl", Dependent  ) ~ "sapling",
    grepl("juven", Dependent  ) ~ "juvenile",
    grepl("mature", Dependent  ) ~ "mature",
    grepl("regen", Dependent  ) ~ "regeneration_pool",
    TRUE ~ "other"
  )) %>%
  mutate(predictor_category = case_when(
    grepl("^manag", Predictor) ~ "Management",  # General SPEI variables
    grepl("^spei", Predictor) ~ "Clim",  # General SPEI variables
    grepl("^drought_spei", Predictor) ~ "Clim",  # Drought-related SPEI
    grepl("^tmp", Predictor) & !grepl("drought", Predictor) ~ "Clim",  # Temperature variables
    grepl("^prcp", Predictor) & !grepl("drought", Predictor) ~ "Clim",  # Precipitation variables
    grepl("^drought_tmp", Predictor) ~ "Clim",  # Drought temperature
    grepl("^drought_prcp", Predictor) ~ "Clim",  # Drought precipitation
    grepl("extract|av\\.nitro", Predictor) ~ "Soil",  # Soil properties like sand, clay, depth, nitrogen
    grepl("distance_edge|disturbance_severity", Predictor) ~ "Disturbance",  # Disturbance variables
    TRUE ~ "Other"  # Catch-all for any other predictor
  )) %>% 
  group_by(Dependent, dependent_category, predictor_category) %>%
  slice_min(AIC, n = 3)   # Select the best 3 based on AIC
# slice(which.max(DevianceExplained))

print(best_predictors)

View(best_predictors)

best_predictors %>% 
  dplyr::filter(dependent_category == 'regeneration_pool')

sjPlot::tab_df(model_metrics,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = FALSE,
               file="outTable/find_best_predictors.doc",
               digits = 1) 


sjPlot::tab_df(best_predictors,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = FALSE,
               file="outTable/best_predictors.doc",
               digits = 1) 




###  variable selection: BORUTA ---------------------------------------------
# List of dependent variables
dependent_vars <- c(#"sum_stems_juvenile", 
                    #"sum_stems_sapling", 
                    "stem_regeneration", # sum of juveniles and saplings
                   # "sum_stems_mature",
                   # "stem_density",      # sum across all classes
                   "advanced",
                   "delayed"
)

# List of predictor variables (spei1 to spei24, drought_spei1 to drought_spei24)
predictor_vars_sub <- c(#"spei1", "spei3", 
                    #"spei6", 
                   # "spei12", 
                    #"spei24",
                    "spei12",
                    "tmp", 
                    #"tmp_z", 
                    "prcp", 
                    #"prcp_z", 
                    #"drought_spei1", "drought_spei3", "drought_spei6", 
                   # "drought_spei12", 
                    #"drought_spei24",
                    #"drought_tmp", 
                    #"drought_prcp",
                    #"salvage_intensity",
                    "management_intensity",
                    # disturbance chars
                    "distance_edge", 
                    "disturbance_severity",
                    
                    # soil
                    #"sand_extract", 
                    "clay_extract", 
                    "depth_extract", 
                    
                    # site info
                    "av.nitro",
                    "richness",
                    'rIVI',
                    "sum_stems_mature",
                    "n_vertical")



## Models: prepare fin tables for individual models  ---------------------------------------------------------------------------------
# test drivers: simplify the analysis:
# Subset the data

df_stem_regeneration2 <- df_fin %>% 
  dplyr::select(all_of(c("stem_regeneration", predictor_vars_sub,"management_intensity",
                                                          "country_pooled", "clim_grid", "clim_class", "x", "y")))


# Centering the variables in your data frame
# Centering all relevant continuous variables in your data frame
df_stem_regeneration2 <- df_stem_regeneration2 %>%
  mutate(
    prcp_c = prcp - mean(prcp, na.rm = TRUE),
    tmp_c = tmp - mean(tmp, na.rm = TRUE),
    spei12_c = spei12 - mean(spei12, na.rm = TRUE),
    distance_edge_c = distance_edge - mean(distance_edge, na.rm = TRUE),
    disturbance_severity_c = disturbance_severity - mean(disturbance_severity, na.rm = TRUE),
    clay_extract_c = clay_extract - mean(clay_extract, na.rm = TRUE),
    av.nitro_c = av.nitro - mean(av.nitro, na.rm = TRUE),
    depth_extract_c = depth_extract - mean(depth_extract, na.rm = TRUE),
    management_intensity_c = management_intensity - mean(management_intensity, na.rm = TRUE)
  )

summary(df_stem_regeneration2)

## Drivers: ---------------------------------------------------------------------


# SIMPLIFY PLOTTING TEST
# Function to generate predictions and plots
create_plot <- function(model, term, data, title, x_label = term, y_label = y_lab, line_color = "blue", fill_color = "blue", scatter = TRUE, x_limit = NULL, scatter_y = "stem_regeneration") {
  
  # Generate predictions
  predicted <- ggpredict(model, terms = term)
  
  # Create the base plot
  plot <- ggplot(predicted, aes(x = x, y = predicted)) +
    geom_line(color = line_color, size = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = fill_color) +
    labs(title = title, x = x_label, y = y_label) +
    theme_classic2()
  
  # Add scatter points if required
  if (scatter) {
    plot <- plot + geom_point(data = data, aes_string(x = term, y = scatter_y), 
                              color = "black", alpha = 0.4, size = 0.5)
  }
  
  # Set x-axis limits if provided
  if (!is.null(x_limit)) {
    plot <- plot + scale_x_continuous(limits = x_limit)
  }
  
  return(plot)
}


# Function to create interaction plots
create_interaction_plot <- function(model, terms, title, data, x_label = terms[1], y_label = y_lab) {
  predicted_interaction <- ggpredict(model, terms = terms)
  
  plot <- ggplot(predicted_interaction, aes(x = x, y = predicted )) +
    geom_line(aes(color = group), size = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
    labs(title = title, x = x_label, y = y_label, color = "Group", fill = "Group") +
    theme_classic2()
  
  return(plot)
}













## Drivers: regeneration density pooled ---------------------------------------------


# Check the structure of the data
str(df_stem_regeneration2)

# Check for missing values
sum(is.na(df_stem_regeneration2))

# Inspect the distribution of each predictor
hist(df_stem_regeneration2$prcp, main="Histogram of Precipitation (prcp)", xlab="prcp")
hist(df_stem_regeneration2$tmp, main="Histogram of Temperature (tmp)", xlab="tmp")

# check corelation between predictors
library(corrplot)

# Select the relevant predictors from your data frame
predictors <- df_stem_regeneration2 %>%
  dplyr::select(prcp_c, tmp_c, spei12_c, distance_edge_c, depth_extract_c, 
                disturbance_severity_c, clay_extract_c, av.nitro_c, management_intensity_c)

# Calculate the correlation matrix
correlation_matrix <- cor(predictors, use = "complete.obs")

# Display the correlation matrix
print(correlation_matrix)
# get basic gam with all preictors
basic_gam <- gam(stem_regeneration     ~ s(prcp_c, k = 5) + s(tmp_c, k = 5) + s(spei12_c, k = 5) + 
                   s(distance_edge_c, k = 5) +
                   s(depth_extract_c, k = 5) +
                   s(disturbance_severity_c, k = 5) +
                   s(clay_extract_c, k = 5) +
                   s(av.nitro_c, k =5) +
                   s(management_intensity, k = 10) +
                   s(country_pooled, bs = 're') # + clim_grid + clim_class
                 , 
                 family = Tweedie(p=1.5), data = df_stem_regeneration2, method = "REML")

int_clim <- gam(stem_regeneration     ~ s(prcp_c, k = 5) + s(tmp_c, k = 5) + s(spei12_c, k = 5) + 
                           s(distance_edge_c, k = 5) +
                           s(depth_extract_c, k = 5) +
                           s(disturbance_severity_c, k = 5) +
                           s(clay_extract_c, k = 5) +
                           s(av.nitro_c, k =5) +
                           s(management_intensity, k = 10) +
                    ti(tmp_c, prcp_c) +
                   s(country_pooled, bs = 're') # + clim_grid + clim_class
                 , 
                 family = Tweedie(p=1.5), data = df_stem_regeneration2, method = "REML")

int_dist <- gam(stem_regeneration     ~ s(prcp_c, k = 5) + s(tmp_c, k = 5) + s(spei12_c, k = 5) + 
                              s(distance_edge_c, k = 5) +
                              s(depth_extract_c, k = 5) +
                              s(disturbance_severity_c, k = 5) +
                              s(clay_extract_c, k = 5) +
                              s(av.nitro_c, k =5) +
                              s(management_intensity, k = 10) +
                             # ti(tmp_c, prcp_c) +
                      ti(disturbance_severity_c, distance_edge_c) +
                      s(country_pooled, bs = 're') # + clim_grid + clim_class
                    , 
                    family = Tweedie(p=1.5), data = df_stem_regeneration2, method = "REML")


AIC(basic_gam,int_clim, int_dist  )
summary(basic_gam)
summary(int_clim)
summary(int_dist)
summary(interaction_model_5_k5)

gam.check(basic_gam)
appraise(basic_gam)
windows()

plot.gam(basic_gam, page = 1)
plot.gam(basic_gam_ti_disturb, page = 1)
plot.gam(basic_gam_ti_clim, page = 1)


# Model with random effect and random slopes for management across countries: account  account for clim_grid to avoid spatial autocorrelation
int_re <- gam(stem_regeneration ~
                             s(prcp_c, k = 5) + s(tmp_c, k = 5) + s(spei12_c, k = 5) + 
                             s(distance_edge_c, k = 5) +
                             s(depth_extract_c, k = 5) +
                             s(disturbance_severity_c, k = 5) +
                             s(clay_extract_c, k = 5) +
                             s(av.nitro_c, k =5) +
                             ti(tmp_c, prcp_c) +
                             ti(disturbance_severity_c, distance_edge_c) +
                             s(management_intensity_c,by = country_pooled, k = 10) + 
                             s(country_pooled, bs = "re") + 
                             s(x,y) +
                             s(clim_grid, bs = "re") 
                             ,
                           family = tw(), method = "REML", data = df_stem_regeneration2)

summary(int_re)
AIC(basic_gam,int_clim, int_dist,int_re ,int_re_dist_spei ,int_re_dist_patch_spei,
    int_re_dist_patch_tmp,
    int_re_dist_prcp,
    int_re_edge_prcp,
    interaction_model_5_k5,
    int_re_dist_prcp_simpl,
    int_re_dist_prcp_simpl2,
    int_re_dist_prcp_simpl3)


# Modify the model to use a random slope for management_intensity by country
int_re_dist_spei <- gam(stem_regeneration ~
                             s(prcp_c, k = 5) + s(tmp_c, k = 5) + 
                             s(spei12_c, k = 5) + 
                             s(distance_edge_c, k = 5) +
                             s(depth_extract_c, k = 5) +
                             s(disturbance_severity_c, k = 5) +
                             s(clay_extract_c, k = 5) +
                             s(av.nitro_c, k =5) +
                             ti(tmp_c, prcp_c, k = 10) +
                             ti(disturbance_severity_c, distance_edge_c, k = 5) +
                             ti(disturbance_severity_c, spei12_c, k = 5) +
                             s(management_intensity_c,by = country_pooled, k = 10) + 
                             s(country_pooled, bs = "re") + 
                             s(x,y) +
                             s(clim_grid, bs = "re") 
                        ,
                           family = tw(), method = "REML", data = df_stem_regeneration2)

# Modify the model to use a random slope for management_intensity by country
int_re_dist_tmp <- gam(stem_regeneration ~
                          s(prcp_c, k = 5) + s(tmp_c, k = 5) + 
                          s(spei12_c, k = 5) + 
                          s(distance_edge_c, k = 5) +
                          s(depth_extract_c, k = 5) +
                          s(disturbance_severity_c, k = 5) +
                          s(clay_extract_c, k = 5) +
                          s(av.nitro_c, k =5) +
                          ti(tmp_c, prcp_c, k = 5) +
                          ti(disturbance_severity_c, distance_edge_c, k = 5) +
                          ti(disturbance_severity_c, tmp_c, k = 5) +
                          s(management_intensity_c,by = country_pooled, k = 10) + 
                          s(country_pooled, bs = "re") + 
                          s(x,y) +
                          s(clim_grid, bs = "re") 
                        ,
                        family = tw(), method = "REML", data = df_stem_regeneration2)


# Modify the model to use a random slope for management_intensity by country
int_re_dist_prcp <- gam(stem_regeneration ~
                         s(prcp_c, k = 5) + s(tmp_c, k = 5) + 
                         s(spei12_c, k = 5) + 
                         s(distance_edge_c, k = 5) +
                         s(depth_extract_c, k = 5) +
                         s(disturbance_severity_c, k = 5) +
                         s(clay_extract_c, k = 5) +
                         s(av.nitro_c, k =5) +
                         ti(tmp_c, prcp_c, k = 5) +
                         ti(disturbance_severity_c, distance_edge_c, k = 5) +
                         ti(disturbance_severity_c, prcp_c, k = 5) +
                         s(management_intensity_c,by = country_pooled, k = 10) + 
                         s(country_pooled, bs = "re") + 
                         s(x,y) +
                         s(clim_grid, bs = "re") 
                       ,
                       family = tw(), method = "REML", data = df_stem_regeneration2)


# simpler: remove the interaction between 

int_re_dist_prcp_simpl <- gam(stem_regeneration ~
                          s(prcp_c, k = 5) + s(tmp_c, k = 5) + 
                          s(spei12_c, k = 5) + 
                          s(distance_edge_c, k = 5) +
                          s(depth_extract_c, k = 5) +
                          s(disturbance_severity_c, k = 5) +
                          s(clay_extract_c, k = 5) +
                          s(av.nitro_c, k =5) +
                          ti(tmp_c, prcp_c, k = 10) +
                         # ti(disturbance_severity_c, distance_edge_c, k = 5) +
                          ti(disturbance_severity_c, prcp_c, k = 10) +
                          s(management_intensity_c,by = country_pooled, k = 10) + 
                          s(country_pooled, bs = "re") + 
                          s(x,y) +
                          s(clim_grid, bs = "re") 
                        ,
                        family = tw(), method = "REML", data = df_stem_regeneration2)



int_re_dist_prcp_simpl2 <- gam(stem_regeneration ~
                                s(prcp_c, k = 5) + s(tmp_c, k = 5) + 
                                s(spei12_c, k = 5) + 
                                s(distance_edge_c, k = 5) +
                                s(depth_extract_c, k = 5) +
                                s(disturbance_severity_c, k = 5) +
                                s(clay_extract_c, k = 5) +
                                s(av.nitro_c, k =5) +
                                ti(tmp_c, prcp_c, k = 10) +
                                 ti(disturbance_severity_c, distance_edge_c, k = 10) +
                                #ti(disturbance_severity_c, prcp_c, k = 10) +
                                s(management_intensity_c,by = country_pooled, k = 4) + 
                                s(country_pooled, bs = "re") + 
                                s(x,y) +
                                s(clim_grid, bs = "re") 
                              ,
                              family = tw(), method = "REML", data = df_stem_regeneration2)


int_re_dist_prcp_simpl3 <- gam(stem_regeneration ~
                                 s(prcp_c, k = 4) + s(tmp_c, k = 4) + 
                                 s(spei12_c, k = 4) + 
                                 s(distance_edge_c, k = 7) +
                                 s(depth_extract_c, k = 4) +
                                 s(disturbance_severity_c, k = 4) +
                                 s(clay_extract_c, k = 5) +
                                 s(av.nitro_c, k =5) +
                                 ti(tmp_c, prcp_c, k = 7) +
                                 #ti(disturbance_severity_c, distance_edge_c, k = 10) +
                                 ti(disturbance_severity_c, prcp_c, k = 7) +
                                 s(management_intensity_c,by = country_pooled, k = 4) + 
                                 s(country_pooled, bs = "re") + 
                                 s(x,y) +
                                 s(clim_grid, bs = "re") 
                               ,
                               family = tw(), method = "REML", data = df_stem_regeneration2)


int_re_dist_prcp_simpl4 <- gam(stem_regeneration ~
                                 s(prcp_c, k = 10) + s(tmp_c, k = 10) + 
                                # s(spei12_c, k = 4) + 
                                 s(distance_edge_c, k = 7) +
                                 #s(depth_extract_c, k = 4) +
                                 s(disturbance_severity_c, k = 4) +
                                 s(clay_extract_c, k = 5) +
                                 #s(av.nitro_c, k =5) +
                                 ti(tmp_c, prcp_c, k = 10) +
                                 #ti(disturbance_severity_c, distance_edge_c, k = 10) +
                                 ti(disturbance_severity_c, prcp_c, k = 5) +
                                 s(management_intensity_c,by = country_pooled, k = 4) + 
                                 s(country_pooled, bs = "re") + 
                                 s(x,y) +
                                 s(clim_grid, bs = "re") 
                               ,
                               family = tw(), method = "REML", data = df_stem_regeneration2)


k.check(int_re_dist_prcp_simpl4)
summary(int_re_dist_prcp_simpl4)

# quick plotting --------------------------------------

df_pred_test <- ggpredict(int_re_dist_prcp_simpl4, 
                          terms = c("disturbance_severity_c", "prcp_c", "country_pooled [AT]"),
                          type = "fixed")

ggplot(df_pred_test, aes(x = x, y = predicted )) +
  geom_line(aes(color = group), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) #+

# add interaction istance_edge and prcp
int_re_edge_prcp <- gam(stem_regeneration ~
                          s(prcp_c, k = 5) + s(tmp_c, k = 5) + 
                          s(spei12_c, k = 5) + 
                          s(distance_edge_c, k = 5) +
                          s(depth_extract_c, k = 5) +
                          s(disturbance_severity_c, k = 5) +
                          s(clay_extract_c, k = 5) +
                          s(av.nitro_c, k =5) +
                          ti(tmp_c, prcp_c, k = 5) +
                          ti(disturbance_severity_c, distance_edge_c, k = 5) +
                          ti(disturbance_severity_c, prcp_c, k = 5) +
                          ti(distance_edge_c, prcp_c, k = 5) +
                          s(management_intensity_c,by = country_pooled, k = 10) + 
                          s(country_pooled, bs = "re") + 
                          s(x,y) +
                          s(clim_grid, bs = "re") 
                        ,
                        family = tw(), method = "REML", data = df_stem_regeneration2)




# interaction istance_edge and spei - not improves moel that much
int_re_dist_patch_spei <- gam(stem_regeneration ~
                          s(prcp_c, k = 5) + s(tmp_c, k = 5) + 
                          s(spei12_c, k = 5) + 
                          s(distance_edge_c, k = 5) +
                          s(depth_extract_c, k = 5) +
                          s(disturbance_severity_c, k = 5) +
                          s(clay_extract_c, k = 5) +
                          s(av.nitro_c, k =5) +
                          ti(tmp_c, prcp_c, k = 5) +
                          ti(disturbance_severity_c, distance_edge_c, k = 5) +
                          ti(disturbance_severity_c, spei12_c, k = 5) +
                            ti(distance_edge_c, spei12_c, k = 5) +
                          s(management_intensity_c,by = country_pooled, k = 10) + 
                          s(country_pooled, bs = "re") + 
                          s(x,y) +
                          s(clim_grid, bs = "re") 
                        ,
                        family = tw(), method = "REML", data = df_stem_regeneration2)

# ad interaction severity with tmp:
int_re_dist_patch_tmp <- gam(stem_regeneration ~
                                s(prcp_c, k = 5) + s(tmp_c, k = 5) + 
                                s(spei12_c, k = 5) + 
                                s(distance_edge_c, k = 5) +
                                s(depth_extract_c, k = 5) +
                                s(disturbance_severity_c, k = 5) +
                                s(clay_extract_c, k = 5) +
                                s(av.nitro_c, k =5) +
                                ti(tmp_c, prcp_c, k = 5) +
                                ti(disturbance_severity_c, distance_edge_c, k = 5) +
                                ti(disturbance_severity_c, spei12_c, k = 5) +
                                ti(distance_edge_c, tmp_c, k = 5) +
                                s(management_intensity_c,by = country_pooled, k = 10) + 
                                s(country_pooled, bs = "re") + 
                                s(x,y) +
                                s(clim_grid, bs = "re") 
                              ,
                              family = tw(), method = "REML", data = df_stem_regeneration2)




summary(int_re_dist_patch_tmp)
interaction_model_5_k5 <- gam(
     stem_regeneration ~ s(prcp, k = 5) + s(tmp, k = 5) + s(spei12, k = 5) +
        s(distance_edge, k = 5) + s(depth_extract, k = 5) + s(disturbance_severity, k = 5) +
         s(clay_extract, k = 5) + s(av.nitro, k = 5) +
         # Change from separate smooths to a random slope model for management_intensity by country
         s(management_intensity, by = country_pooled, bs = "re") +
         s(country_pooled, bs = "re") +  # Random intercept for country
         ti(prcp, tmp, k = 10) + 
        #ti(spei12, tmp, k = 10) +
        s(x, y) + 
        s(clim_grid, bs = "re"), 
    family = tw(), 
     method = "REML", 
    data = df_stem_regeneration2)


interaction_model_5_k5_no_spei <- gam(
  stem_regeneration ~ s(prcp, k = 5) + s(tmp, k = 5) +# s(spei12, k = 5) +
    s(distance_edge, k = 5) + s(depth_extract, k = 5) + s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) + s(av.nitro, k = 5) +
    # Change from separate smooths to a random slope model for management_intensity by country
    s(management_intensity, by = country_pooled, bs = "re") +
    s(country_pooled, bs = "re") +  # Random intercept for country
    ti(prcp, tmp, k = 10) + 
    #ti(spei12, tmp, k = 10) +
    s(x, y) + 
    s(clim_grid, bs = "re"), 
  family = tw(), 
  method = "REML", 
  data = df_stem_regeneration2)



# store the best model for regeneration density
fin.m.reg.density <- int_re_dist_prcp_simpl3 

vis.gam(fin.m.reg.density, view = c("prcp_c", "tmp_c"), plot.type = "persp",
        main = "Interaction between Precipitation and Temperature",
        zlab = "Stem Regeneration", xlab = "Precipitation", ylab = "Temperature")

# ti(disturbance_severity, distance_edge): Interaction between disturbance severity and distance to edge
vis.gam(fin.m.reg.density, view = c("disturbance_severity_c", "distance_edge_c"), plot.type = "persp",
        zlab = "Stem Regeneration")

# ti(disturbance_severity, distance_edge): Interaction between disturbance severity and distance to edge
vis.gam(fin.m.reg.density, view = c("disturbance_severity_c", "prcp_c"), plot.type = "persp",
        zlab = "Stem Regeneration")


# check again correlation between predictors:

# Load necessary libraries
library(corrplot)

# Select predictors from your data frame
predictors <- df_stem_regeneration2[, c("prcp_c", "tmp_c", "spei12_c", 
                                        "distance_edge_c", "depth_extract_c",
                                        "disturbance_severity_c", 
                                        "clay_extract_c", "av.nitro_c", 
                                        "management_intensity_c")]

# Calculate correlation matrix
correlation_matrix <- cor(predictors,  method = "spearman", use = "complete.obs")

# Plot the correlation matrix
corrplot(correlation_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, 
         title = "Correlation Matrix of Predictors",
         addCoef.col = "black", number.cex = 0.7)


#summary(interaction_model_4)
appraise(fin.m.reg.density)


# test for spatial autocorrelation

# Extract model residuals
model_residuals <- residuals(fin.m.reg.density, type = "pearson")

# Load necessary libraries
library(spdep)

# Create coordinates matrix
coords <- cbind(df_stem_regeneration2$x, df_stem_regeneration2$y)

# Create a spatial neighbors object (e.g., using k-nearest neighbors)
# Adjust k based on the density and distribution of your data points
nb <- knn2nb(knearneigh(coords, k = 5))

# Convert neighbors list to a weights list
listw <- nb2listw(nb, style = "W")
# Perform Moran's I test
moran_test <- moran.test(model_residuals, listw)
print(moran_test)


# PLot  ---------------------------------------------------------------------------
y_lab = 'Stem density [#/ha]'

summary(fin.m.reg.density)

# Create plots for individual variables
plot_prcp <- create_plot(fin.m.reg.density, "prcp_c", df_stem_regeneration2, "prcp**", line_color = "blue", fill_color = "blue")
plot_tmp <- create_plot(fin.m.reg.density, "tmp_c", df_stem_regeneration2, "tmp**", line_color = "red", fill_color = "red")
plot_disturbance_severity <- create_plot(fin.m.reg.density, "disturbance_severity_c", df_stem_regeneration2, "Dist Severity*", line_color = "purple", fill_color = "purple")
plot_distance_edge <- create_plot(fin.m.reg.density, "distance_edge_c", df_stem_regeneration2, "Dist edge 0.12", x_limit = c(0, 600), line_color = "grey", fill_color = "grey")
(plot_disturbance_severity)

# Create interaction plots
plot_interaction1 <- create_interaction_plot(fin.m.reg.density, 
                                             c("prcp_c", "tmp_c"), "prcp & tmp 0.05", df_stem_regeneration2) +
  labs(color = "tmp", fill = "tmp") 
(plot_interaction1)

plot_interaction2 <- create_interaction_plot(fin.m.reg.density, 
                                             c("disturbance_severity_c", "distance_edge_c"), "lalal", df_stem_regeneration2) +
  labs(color = "disturbance_severity_c", fill = "distance_edge_c") 
(plot_interaction2)

plot_interaction3 <- create_interaction_plot(fin.m.reg.density, 
                                             c("disturbance_severity_c", "prcp_c"), "lala", df_stem_regeneration2) +
  labs(color = "disturbance_severity_c", fill = "prcp_c") 
(plot_interaction3)


# uderstand why data are different

predicted_interaction <- ggpredict(fin.m.reg.density, 
                                   terms = c("disturbance_severity_c", "prcp_c"),
                                   type = "fixed")

ggplot(predicted_interaction, aes(x = x, y = predicted )) +
  geom_line(aes(color = group), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) #+
  #' labs(title = title, x = x_label, y = y_label, color = "Group", fill = "Group") +
  #' theme_classic2()

# Combine the individual plots into one figure




library(cowplot)
#plot_grid(p1, p2, labels = "auto", label_size = 12)
comb_upper  <- plot_grid(plot_interaction1, plot_interaction2,ncol = 2, labels = c("[a]","[b]"))
comb_middle <- plot_grid(plot_mgmt_intensity, plot_countries,ncol = 2,labels = c("[c]","[d]"))
comb_lower <- plot_grid(plot_prcp, plot_tmp,
                            plot_spei12, plot_disturbance_severity, plot_distance_edge, 
                            plot_management_intensity, ncol = 3,
                        labels = c("[e]","[f]","[g]","[h]","[i]","[j]"))
# Combine the individual plots into one figure
combined_plot <- plot_grid(comb_upper,
                           comb_middle, 
                           comb_lower,
                           ncol = 1, nrow = 3,
                           rel_heights = c(1.5, 2,3))

combined_plot


# Save the combined plot
ggsave('outFigs/fig_regen_pool_drivers.png', plot = combined_plot, width = 7, height = 9.5, bg = 'white')



# Identify random effects using the model's "smooth" component
smooth_terms <- summary(fin.m.advanced)$s.table

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
sjPlot::tab_model(fin.m.advanced,
                  show.re.var = TRUE,        # Show the variance components
                  #show.icc = TRUE,           # Show Intraclass Correlation Coefficient
                  #show.dev = TRUE,           # Show deviance
                  pred.labels = c("Intercept", pred_labels), # Replace smooth term labels
                  dv.labels = paste0("Explained Deviance: ", round(100 * summary(model)$dev.expl, 2), "%"), 
                  file = "outTable/full_drivers_advanced.doc")











## Wilcox: Boxplot sites differences: delayed vs advaced ----------------------------

# two categories: count how many plots i have?

prop.table(table(df_fin$adv_delayed_wider))
prop.table(table(df_fin$adv_delayed))


# Select relevant columns including 'RegenerationStatus' and the desired variables
variables_to_plot <- c(
  "tmp",
  "prcp",     
  "spei12", 
  "stem_regeneration",
  "sum_stems_juvenile"  ,
  "sum_stems_sapling", 
  "sum_stems_mature",
  "av.nitro", 
  "depth_extract", 
  "management_intensity", 
  "salvage_intensity", 
  "protection_intensity", 
  "distance_edge" , 
  "disturbance_severity",
  "salvage_intensity",
  "clay_extract"#,
 # "adv_delayed"
                       )

# differentiate classes:

# Step 2: Calculate median and IQR for each variable by regeneration status
summary_stats_narrow <- 
  df_fin %>%
  na.omit() %>% 
    dplyr::select(all_of(c(variables_to_plot, 'adv_delayed'))) %>% 
  gather(key = "Variable", value = "Value", -adv_delayed) %>%
  group_by(adv_delayed, Variable) %>%
  summarise(
    Median = median(Value, na.rm = TRUE),
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE)
  ) %>%
  mutate(IQR_Lower = Median - Q1, 
         IQR_Upper = Q3 - Median)


# make wider categories: for delayed and advanced

summary_stats_wider <- 
  df_fin %>%
  na.omit() %>% 
    dplyr::select(all_of(c(variables_to_plot, 'adv_delayed_wider'))) %>% 
  gather(key = "Variable", value = "Value", -adv_delayed_wider) %>%
  group_by(adv_delayed_wider, Variable) %>%
  summarise(
    Median = median(Value, na.rm = TRUE),
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE)
  ) %>%
  mutate(IQR_Lower = Median - Q1, 
         IQR_Upper = Q3 - Median)



# chack my categories?
df_fin %>% 
  dplyr::select(site, 
                stem_density,
                stem_regeneration,
                sum_stems_juvenile,
                 sum_stems_sapling, 
                 sum_stems_mature,
                adv_delayed,
                adv_delayed_wider
                ) %>% 
  View()


# Step 3: Create the median and IQR plot using ggplot2
ggplot(summary_stats_narrow, aes(x = Variable, y = Median, color = adv_delayed)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Points for median
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, 
                position = position_dodge(width = 0.5)) +          # Error bars for IQR
  labs(title = "Median and IQR Plot for Delayed vs Advanced Regeneration",
       x = "Variables", y = "Median and IQR") +
  theme_classic() +
  facet_wrap(.~Variable, scales = 'free')+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_color_manual(values = c("blue", "red")) +  # Set colors for delayed vs advanced
  theme(legend.title = element_blank())            # Remove legend title

# for wider interavls:

# Step 3: Create the median and IQR plot using ggplot2
ggplot(summary_stats_wider , aes(x = Variable, y = Median, color = adv_delayed_wider  )) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Points for median
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, 
                position = position_dodge(width = 0.5)) +          # Error bars for IQR
  labs(title = "Median and IQR Plot for Delayed vs Advanced Regeneration: wider intervals",
       x = "Variables", y = "Median and IQR") +
  theme_classic() +
  facet_wrap(.~Variable, scales = 'free')+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_color_manual(values = c("blue", "red")) +  # Set colors for delayed vs advanced
  theme(legend.title = element_blank())            # Remove legend title




### test differences between groups: Narrow --------------------------------------
# Step 1: Reshape the data to long format
df_long_narrow <- df_fin %>%
  na.omit() %>% 
  dplyr::select(tmp, 
                prcp, 
                spei12,
                disturbance_severity, 
                distance_edge, 
                adv_delayed) %>% 
  #dplyr::select(all_of(variables_to_plot), adv_delayed) %>% 
  gather(key = "Variable", value = "Value", -adv_delayed)


# list groups to pairwise comparison
comparisons <- list(c("Delayed", "Other"), c("Delayed", "Advanced"), c("Other", "Advanced"))

# Plot using ggboxplot
p_boxplot_wilcox_narrow <- ggboxplot(df_long_narrow, x = "adv_delayed", y = "Value", 
                              fill = "adv_delayed", 
                              palette = c("blue", "red", "green"),
                              facet.by = "Variable", scales = "free_y", 
                              ylab = "Values", xlab = "Regeneration Status",
                              outlier.size = .2,
                              size = 0.2) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", 
                     #label.y = label_y_positions[as.character(df_long$Variable)], # Use the calculated y positions
                     size = 2,
                     label.x = 1.5) +  # Position labels between the groups+
  labs(title = "Delayed <=50, advanced >= 1000 juv",
       x = "Reg. Status", y = "Vals")+
  theme(
    legend.position = 'none',
    text = element_text(size = 3),         # Set all text to size 3
    axis.text = element_text(size = 3),    # Axis tick labels
    axis.title = element_text(size = 3),   # Axis titles
    strip.text = element_text(size = 3),   # Facet labels
    legend.text = element_text(size = 3),  # Legend text
    plot.title = element_text(size = 5),   # Plot title
    strip.background = element_blank(),    # Remove the box around facet names
    strip.placement = "outside",           # Optional: Move facet label outside the plot area
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)  # Add a square border around the plots
  ) 

p_boxplot_wilcox_narrow

# Save the combined plot (optional)
# Save the plot ensuring text sizes are preserved
ggsave("outFigs/p_boxplot_wilcox_delayed_50.png", plot = p_boxplot_wilcox_narrow, 
       width = 3, height = 3.2, units = "in", dpi = 300, 
       bg = 'white', scale = 1)

### wider intervals: ------------------------------------------------------

# Step 1: Reshape the data to long format
df_long_wider <- df_fin %>%
  dplyr::select(tmp, 
                prcp, 
                spei12,
                disturbance_severity, 
                distance_edge, 
                adv_delayed_wider) %>% 
  na.omit() %>% 
 # dplyr::select(all_of(variables_to_plot), adv_delayed_wider) %>% 
  gather(key = "Variable", value = "Value", -adv_delayed_wider)



# get pairwised comparison:


p_boxplot_wilcox_wider <- ggboxplot(df_long_wider, x = "adv_delayed_wider", y = "Value", 
                              #color = "adv_delayed_wider",
                              fill = "adv_delayed_wider",
                              palette = c("blue", "red", "green"),
                              facet.by = "Variable", scales = "free_y", 
                              outlier.size = .2,
                              size = 0.2,
                              ylab = "Values", xlab = "Regeneration Status") +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", 
                     #label.y = label_y_positions[as.character(df_long$Variable)], # Use the calculated y positions
                     size = 2,
                     label.x = 1.5) +  # Position labels between the groups
  #geom_boxplot(lwd = 0.3) +      
  labs(title = "Delayed <=500, advanced >= 1000 juv",
       x = "Reg. Status", y = "Vals") +
  theme(
    legend.position = 'none',
    text = element_text(size = 3),         # Set all text to size 3
    axis.text = element_text(size = 3),    # Axis tick labels
    axis.title = element_text(size = 3),   # Axis titles
    strip.text = element_text(size = 3),   # Facet labels
    legend.text = element_text(size = 3),  # Legend text
    plot.title = element_text(size = 5),   # Plot title
    strip.background = element_blank(),    # Remove the box around facet names
    strip.placement = "outside",           # Optional: Move facet label outside the plot area
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)  # Add a square border around the plots
  ) 


p_boxplot_wilcox_wider

# Save the combined plot (optional)
ggsave("outFigs/p_boxplot_wilcox_delayed_500.png", p_boxplot_wilcox_wider, width = 3, height = 3.2, 
       dpi = 300,  bg = 'white')


#### Wilcox tablle -----------------------------------------------
library(tidyr)

# Perform pairwise Wilcoxon test
pairwise_results <- df_long_wider %>%
  group_by(Variable) %>%
  summarise(pairwise_p = list(pairwise.wilcox.test(Value, adv_delayed_wider, p.adjust.method = "BH")$p.value))

# Expand pairwise Wilcoxon p-value matrix into a long format
pairwise_table <- pairwise_results %>%
  mutate(pairwise_p = map(pairwise_p, ~as.data.frame(as.table(.)))) %>% # Convert matrix to data frame
  unnest(pairwise_p) %>%
  rename(Comparison1 = Var1, Comparison2 = Var2, p_value = Freq) %>%
  mutate(Comparison = paste(Comparison1, "vs", Comparison2)) %>%
  dplyr::select(Variable, Comparison, p_value)

# Show the pairwise Wilcoxon p-values as a table
pairwise_table


# Assuming `df_delayed_filtered` is your data frame
# Splitting data into delayed and advanced groups based on stem_regeneration
df_delayed <- df_fin[df_fin$adv_delayed == "Delayed", ]
df_advanced <- df_fin[df_fin$adv_delayed == "Advanced", ]

# List of continuous variables to test
variables <- c("rIVI", "sum_stems_mature", "prcp", "spei12", "drought_spei12", "spei3", "drought_spei3", 
               "tmp", "av.nitro", "depth_extract", "management_intensity")

# Apply t-test or Mann-Whitney U test
test_results <- lapply(variables, function(var) {
  delayed_values <- df_delayed[[var]]
  advanced_values <- df_advanced[[var]]
  
  # Check if data is normally distributed
  normality_delayed <- shapiro.test(delayed_values)$p.value
  normality_advanced <- shapiro.test(advanced_values)$p.value
  
  if (normality_delayed > 0.05 & normality_advanced > 0.05) {
    test <- t.test(delayed_values, advanced_values)  # Use t-test if data is normal
  } else {
    test <- wilcox.test(delayed_values, advanced_values)  # Use Mann-Whitney U test if not normal
  }
  
  return(data.frame(
    Variable = var,
    Test = ifelse(normality_delayed > 0.05 & normality_advanced > 0.05, "T-test", "Mann-Whitney"),
    p_value = test$p.value
  ))
})

# Combine results into a single data frame
test_results_df <- do.call(rbind, test_results)
test_results_df$p_value <- round(test_results_df$p_value,2 )
print(test_results_df)


# Removing rows with missing data (if any) -------------------------------------
df_delayed <- na.omit(df_delayed)
df_advanced <- na.omit(df_advanced)
df_stem_regeneration <- na.omit(df_stem_regeneration)



# my data are of course aggregated on climatic lelevl!!! (clim vars from ERA NET!)

# average dependent variable bassed on unique tmp and prec
# Step 1: Create a unique ID for each unique combination of tmp and prcp

# Step 2: Aggregate stem_density and other necessary columns by this unique ID
df_fin_agg <- df_fin %>%
  group_by(clim_id) %>%
  summarise(
    avg_stem_regeneration = mean(stem_regeneration, na.rm = TRUE),  # Average stem density
    tmp_c = first(tmp_c),  # Retain the first tmp value for each group (all are the same in the group)
    prcp_c = first(prcp_c),  # Retain the first prcp value
    tmp = first(tmp),  # Retain the first tmp value for each group (all are the same in the group)
    prcp = first(prcp),  # Retain the first prcp value
    richness = mean(richness, na.rm = TRUE),  # Example: average other columns if needed
    management_intensity = mean(management_intensity, na.rm = TRUE)  # Average other variables if required
  ) %>%
  ungroup()


df_fin_agg <- df_fin %>%
  group_by(clim_id) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE),  # Average all numeric columns
            across(where(is.character), first),             # Retain the first value for all character columns
            across(where(is.factor), first))   


# Fit your model with the aggregated data
interaction_model_agg <- gam(stem_regeneration ~ 
                               s(tmp_c) + 
                               s(prcp_c) + 
                               ti(tmp_c, prcp_c) + 
                               # Include other variables as needed
                               s(richness) + 
                               s(management_intensity),
                             family = Tweedie(p = 1.46), 
                             method = 'REML',
                             data = df_fin_agg)

summary(interaction_model_agg)
appraise(interaction_model_agg)

hist(df_fin_agg$stem_regeneration)


# Visualize the main effect of tmp_c
library( ggeffects)
tmp_effect <- ggpredict(interaction_model_smooth, terms = "tmp_c")
p1 <- plot(tmp_effect) + ggtitle("Effect of Temperature on Stem Regeneration")

# Visualize the main effect of prcp_c
prcp_effect <- ggpredict(interaction_model_smooth, terms = "prcp_c")
p2 <- plot(prcp_effect) + ggtitle("Effect of Precipitation on Stem Regeneration")

# Visualize the interaction effect of tmp_c and prcp_c
interaction_effect <- ggpredict(interaction_model_smooth, terms = c("tmp_c", "prcp_c"))
p3 <- plot(interaction_effect) + ggtitle("Interaction: Temperature and Precipitation on Stem Regeneration")
ggarrange(p1,p2,p3)



# understand smooths vs ggpredict
m_c <- gam(stem_regeneration ~ #ti(x, y) + 
      s(tmp_c) +              # Smooth for temperature
      s(prcp_c) +             # Smooth for precipitation
      ti(tmp_c, prcp_c),# +     # Smooth interaction of temperature and precipitation
      family = Tweedie(p = 1.46), 
    method = 'REML',
    data = df_fin)


m <- gam(stem_regeneration ~ #ti(x, y) + 
             s(tmp) +              # Smooth for temperature
             s(prcp) +             # Smooth for precipitation
             ti(tmp, prcp),# +     # Smooth interaction of temperature and precipitation
           family = Tweedie(p = 1.46), 
           method = 'REML',
           data = df_fin)


AIC(m_c, m)
draw(m_c)
draw(m)


# plotting :
tmp_effect <- ggpredict(m_c, terms = "tmp_c")
p1_c <- plot(tmp_effect) + ggtitle("Effect of Temperature on Stem Regeneration")

# Visualize the main effect of prcp_c
prcp_effect <- ggpredict(m_c, terms = "prcp_c")
p2_c <- plot(prcp_effect) + ggtitle("Effect of Precipitation on Stem Regeneration")

# Visualize the interaction effect of tmp_c and prcp_c
interaction_effect <- ggpredict(m_c, terms = c("tmp_c", "prcp_c"))
p3_c <- plot(interaction_effect) + ggtitle("Interaction: Temperature and Precipitation on Stem Regeneration")
ggarrange(p1_c,p2_c,p3_c)


# no centered, simple: --------------------- 
# plotting :
tmp_effect <- ggpredict(m, terms = "tmp")
p1 <- plot(tmp_effect) + ggtitle("Effect of Temperature on Stem Regeneration")

# Visualize the main effect of prcp_c
prcp_effect <- ggpredict(m, terms = "prcp")
p2 <- plot(prcp_effect) + ggtitle("Effect of Precipitation on Stem Regeneration")

# Visualize the interaction effect of tmp_c and prcp_c
interaction_effect <- ggpredict(m, terms = c("tmp", "prcp"))
p3 <- plot(interaction_effect) + ggtitle("Interaction: Temperature and Precipitation on Stem Regeneration")
ggarrange(p1,p2,p3)









# Step 9: Add random effects for region and country
random_effects_model <- gam(stem_regeneration ~ ti(x, y) + s(tmp, k = 10) + s(prcp, k = 10) + 
                              s(distance_edge, k = 10) + s(disturbance_severity) + 
                              s(management_intensity, by = country_pooled, k = 5) + 
                              s(av.nitro, k = 15) + s(clay_extract, k = 10) + 
                              ti(tmp, prcp) + s(region, bs = 're') + s(country_abbr, bs = 're'), 
                            family = Tweedie(p = 1.46), 
                            method = 'REML',
                            data = df_fin)
summary(random_effects_model)



# Step 9: Add random effects for region and country
random_effects_model_simpl <- gam(stem_regeneration ~ ti(x, y) + s(tmp, k = 10) + s(prcp, k = 10) + 
                              s(distance_edge, k = 10) + s(disturbance_severity) + 
                              s(management_intensity, k = 5) +  # , by = country_pooled, 
                              s(av.nitro, k = 15) + s(clay_extract, k = 10) + 
                              ti(tmp, prcp) + s(region, bs = 're') + s(country_abbr, bs = 're'), 
                            family = Tweedie(p = 1.46), 
                            method = 'REML',
                            data = df_fin)
summary(random_effects_model_simpl)




# Step 10: Compare models using AIC
AIC_comparison <- AIC(base_model, drought_model, temp_model, precip_model, distance_model, disturbance_model, 
                      management_model, soil_model, interaction_model, random_effects_model, random_effects_model_simpl)
print(AIC_comparison)

# Step 11: Check residual plots and diagnostics for the final model
best_model <- interaction_model  # Choose the best model based on AIC and performance

# Visualize residuals and fitted values
appraise(best_model)
plot(best_model, page = 1)

# Step 12: Check for concurvity (multicollinearity)
concurvity(best_model)

# Step 13: Perform k-check to validate smoothers (overfitting)
k.check(best_model)

# Final model diagnostics:
summary(best_model)










# Model selection ---------------------------------------------------------------
library(MuMIn)

### drivers juveniles -----------------------------------------------------------------------------
# Fit a global model with all possible predictors
global_model_juv <- gam(sum_stems_juvenile ~ #s(drought_spei6, k = 8) +
                          s(tmp_z, k = 7) +
                          s(prcp_z, k = 7) +
                          s(drought_spei12, k = 10) +
                          #s(drought_spei24, k = 10) +
                          s(distance_edge, k = 10) +  # did not improved the mode
                          disturbance_severity +
                          clay_extract  +
                          s(clay_extract, k = 10) +
                          s(depth_extract, k = 10) +
                          s(av.nitro, k = 15) +
                          s(country_abbr, bs = "re"),
                        family = tw(), 
                        data = df_fin,
                        method = "REML",
                        na.action = na.fail)

concurvity(global_model_juv)


# the negative binomail performs way less then tw distribution.
model_nb <- gam(sum_stems_juvenile ~  s(tmp_z, k = 7) +
                  s(prcp_z, k = 7) +
                  s(drought_spei12, k = 10) +
                  #s(drought_spei24, k = 10) +
                  s(distance_edge, k = 10) +  # did not improved the mode
                  disturbance_severity +
                  clay_extract  +
                  s(clay_extract, k = 10) +
                  s(depth_extract, k = 10) +
                  s(av.nitro, k = 15) +
                  s(country_abbr, bs = "re"),
                family = nb(), data = df_fin)

####  try to adjust tw parameters -------------------------------------------

# optimize tweedie values -------------------------------------------

# Define a grid of p values that are more appropriate for data with zeros
p_values <- seq(1.4, 1.6, by = 0.01)  # Focus on p-values between 1 and 2 to account for presence to 0
# in general p = 1 - closer to poisson (can be zero inflated) 
#            p = 2 closer to gamma (can only have positive values)

# Initialize a placeholder to store results
results <- list()

# Loop over different p values
for (p_val in p_values) {
  model <- gam(sum_stems_juvenile ~ s(tmp_z, k = 7) + s(prcp_z, k = 7) + 
                 #s(drought_spei12, k = 10) +
                 s(distance_edge, k = 10) +
                 disturbance_severity + clay_extract + s(clay_extract, k = 10) +
                 s(depth_extract, k = 10) + s(av.nitro, k = 15) +
                 s(country_abbr, bs = "re") + 
                 ti(x,y), 
               family = Tweedie(p = p_val),  # Use p between 1 and 2
               data = df_fin)
  
  # Store model AIC and the value of p
  results[[as.character(p_val)]] <- AIC(model)
}

# Find the best p value (lowest AIC)
best_p <- as.numeric(names(which.min(unlist(results))))
print(paste("Optimal p value:", best_p))

# best is p = 1.5 orp =  1.46














# vertical structure vs disturbance patch size? ---------------------------

df_fin %>% 
  ggplot(aes(x = distance_edge,
             y  = n_vertical)) +
  geom_jitter() +
  geom_smooth() +
  theme_classic()



# Country effect ----------------------------------------------------------


df_richness <- df_richness %>% 
  dplyr::filter(cluster %in% cluster_id_keep) 


p.country.richness <- ggplot(data = df_richness,
                             aes(x = country, y = richness,fill = country)) + 
  # Add bars as medians
  stat_summary(fun = "median", 
               geom = "bar") +
  stat_summary(
    data = df_richness,
    mapping = aes(x = country, y = richness),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) +
  ylab('Richness [#]') + 
  theme_classic() +
  theme(legend.position = 'none')
#coord_cartesian(ylim=c(60, 67)) # ylim=c(59,66)



stem_dens_ha_cluster_sum <- stem_dens_ha_cluster_sum %>% 
  dplyr::filter(cluster %in% cluster_id_keep) 


df_stems_country <- stem_dens_ha_cluster_sum %>% 
  ungroup(.) %>% 
  group_by(cluster, manag, country) %>% # sum stems per luster, manag and country
  dplyr::reframe(sum_stems = sum(total_stems)) 

p.country.density <-  df_stems_country %>% 
  ggplot(aes(x = country, y = sum_stems,fill = country )) +
  # Add bars as medians
  stat_summary(fun = "median", 
               geom = "bar") +
  stat_summary(
    data = df_stems_country,
    mapping = aes(x = country, y = sum_stems),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) +
  ylab('Stems [ha]') + 
  theme_classic() +
  theme(legend.position = 'none')

(p.country.density)
#coord_cartesian(ylim=c(60, 67)) # ylim=c(59,66)

# test significance - for 2 groups only
wilcox.test(sum_stems ~ manag, data = df_stems)


# for several groups:
kruskal.test(sum_stems ~ country, data = df_stems_country)


# Perform the Dunn test
dunnTest <- dunn.test(df_stems_country$sum_stems, df_stems_country$manag, method="bonferroni")
dunnTest


ggarrange(p.country.density, p.country.richness)


# Save models --------------------------------------------------

# Save the model object and input data
save(#fin.m.delayed, df_delayed2, 
     fin.m.advanced, df_advanced2, 
     
     #fin.m.reg.density, df_stem_regeneration2,
     file = "outData/stem_density_models.RData")

# 4. future developmenet -------------------------------------------

# compare current species composition with future ones form wessely

# 
head(stem_dens_species_long_cluster)

# get species merging table
species_look_up_simple<- read.csv("rawData/tree_sp_simple.csv", sep = ';')

# get species merging table
species_look_up_full<- read.csv("rawData/tree_sp_field_wessely_merged.csv", sep = ';')


# identify what species are present per cluster
present_species <- 
  stem_dens_species_long_cluster %>% 
  ungroup() %>% 
  group_by(cluster, Species) %>% 
  summarize(sum_stem_density = sum(stem_density, na.rm = T )) %>% 
  mutate(presence = if_else(sum_stem_density > 0, 1, 0)) %>% 
  dplyr::select(-sum_stem_density ) %>% 
  left_join(species_look_up_simple, by = c("Species" = "acc")) %>% 
  dplyr::select(-latin) 

present_species <- present_species %>% 
  dplyr::rename(acc = Species) %>% 
  dplyr::rename(site = cluster)

# read table with Wessely species
future_species <- fread('outTable/species_presence_clim_change.csv')
future_species_full <- fread('outTable/species_presence_clim_change_full.csv')




# add acronyms and consider presence across several species: eg betula sp.
future_species_sum <- 
  future_species %>%
  left_join(species_look_up_full, by = c('species' = 'wessely')) %>%
  # group by species to allow occurence of species that have specified genus: eg betula sp.
  group_by(site, scenario, acc) %>% 
  # Summarize and set sum_presence to 1 if the sum is greater than 1
  summarize(sum_presence = pmin(sum(overall_presence), 1), .groups = 'drop')
  
wide_future_species <- future_species_sum %>%
  pivot_wider(names_from = scenario, values_from = sum_presence) 

# merge both tables: the presently recorded species and species under climate scenarios
df_compare_future_species <- wide_future_species %>% 
  full_join(present_species) %>% 
  dplyr::rename(current = presence)

# add country indication
df_compare_future_species <- df_compare_future_species %>%
  # Extract the first two characters of 'site' as 'region' and convert to integer
  mutate(region = as.integer(substr(site, 1, 2))) %>%
  # Left join with unique_regions_per_country to get country indication
  left_join(unique_regions_per_country, by = c("region" = "unique_regions"))


# how many plots per cpuountry does not contain any currently present species??? ---------------
# evaluate species on plot level: Share of plots where none of the currently present species 
# remain within their climate niche until the end of the 21st century

# check for single country: how many plots are there without any current species present?

# seems correct
  df_compare_future_species %>%
  dplyr::filter(country_pooled == "SK") %>% 
  group_by(site, country_pooled) %>%
    dplyr::filter(site == "16_115") %>% 
    View()
  # summarise(
  #   # Species that remain present (presence == 1 and also present in any future scenario)
  #   same_rcp26 = sum(current == 1 & rcp26 == 1, na.rm = T ), # some nAs can be preent if species was in the field (pseudotsusa mensioes : psme, but does not exists in Wesely database)
  #   same_rcp45 = sum(current == 1 & rcp45 == 1, na.rm = T ),
  #   same_rcp85 = sum(current == 1 & rcp85 == 1, na.rm = T ) #
  # ) %>% 
  # ungroup() 






# Compare presence with future scenarios: caount the number of plots that have not a single species spared between curremt and future scenarios
df_compare_future_species_summary <-
  
  df_compare_future_species %>%
  group_by(site, country_pooled) %>%
  summarise(
    # Species that remain present (presence == 1 and also present in any future scenario)
    same_rcp26 = sum(current == 1 & rcp26 == 1, na.rm = T ), # some nAs can be preent if species was in the field (pseudotsusa mensioes : psme, but does not exists in Wesely database)
    same_rcp45 = sum(current == 1 & rcp45 == 1, na.rm = T ),
    same_rcp85 = sum(current == 1 & rcp85 == 1, na.rm = T ) #
  ) %>% 
    ungroup() 

# Print the result
head(df_compare_future_species_summary)


plot_summary_share <- 
  df_compare_future_species_summary %>%
  group_by(country_pooled) %>%
  summarize(
    total_plots = n(),
    plots_none_present_rcp26 = sum(same_rcp26 < 1, na.rm = TRUE),
    plots_none_present_rcp45 = sum(same_rcp45 < 1, na.rm = TRUE),
    plots_none_present_rcp85 = sum(same_rcp85 < 1, na.rm = TRUE)
  ) %>% 
  mutate(share_none_present26 = plots_none_present_rcp26/total_plots*100,
         share_none_present45 = plots_none_present_rcp45/total_plots*100,
         share_none_present85 = plots_none_present_rcp85/total_plots*100)

# View the result
print(plot_summary_share)

# Calculate shares (percentages) and format output for each RCP scenario
plot_summary_formatted_share <- plot_summary_share %>%
  mutate(
    rcp26 = paste0(plots_none_present_rcp26, " (", round((plots_none_present_rcp26 / total_plots) * 100, 1), "%)"),
    rcp45 = paste0(plots_none_present_rcp45, " (", round((plots_none_present_rcp45 / total_plots) * 100, 1), "%)"),
    rcp85 = paste0(plots_none_present_rcp85, " (", round((plots_none_present_rcp85 / total_plots) * 100, 1), "%)")
  ) %>%
  dplyr::select(country_pooled, total_plots, rcp26, rcp45, rcp85)


print(plot_summary_formatted_share)

# Evaluate by species: what species are present/ country? and how many of them are not clim suitbale?

# Test for single country: 


current_species<- 
  df_compare_future_species %>%  
  dplyr::filter(country_pooled == "SK") %>% 
  dplyr::filter(current == 1) %>% 
  #dplyr::filter(rcp26 == 1) %>% 
  ungroup() %>% 
  distinct(acc) #%>%
 # rename(current = acc)
  #pull()
 
rcp26_species<- df_compare_future_species %>%  
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
  df_compare_future_species %>%  
  dplyr::filter(country_pooled == "SK") %>% 
  dplyr::filter(current == 1) %>% 
  ungroup() %>% 
  distinct(acc) %>%
 pull()

rcp26_species_v<- df_compare_future_species %>%  
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


4/12
x = c('a', 'b', 'c')
y = c('a','c', 'd', 'e', 'f')

sp_diff <- setdiff(x, y)
length(sp_diff)/length(x)*100



# run for all species: -----------------------------------------------------



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
species_loss_summary <- df_compare_future_species %>%
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
