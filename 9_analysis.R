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



# Read data -----------------------------------------------------------------------

# get vegetation data
load("outData/veg.Rdata")

# final tables on site level
df_fin <- fread('outData/indicators_for_cluster_analysis.csv')

# ad disturbance severity based on presence/absence of mature trees 
df_mature_dist_severity <- fread('outData/disturb_severity_mature.csv')

df_fin <- df_fin %>% 
  left_join(df_mature_dist_severity, by = c('site' = 'cluster'))

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
   mutate(stem_regeneration = sum_stems_juvenile + sum_stems_sapling)


# export as xy coordinates climate cluster categorized 
df_fin_clim_clust_xy <- st_as_sf(df_fin, coords = c("x", "y"), crs = crs(xy))  

# Step 2: Check the structure of the sf object to ensure everything is correct
print(df_fin_clim_clust_xy)

# Step 3: Export the data as a GeoPackage (GPKG)
#st_write(df_fin_clim_clust_xy, "outData/xy_clim_cluster.gpkg", layer = "df_fin", driver = "GPKG", append = FALSE)

# get only a ataframe of teh climate clusters and sites for easy merging to detiailed veg data
#clim_cluster_indicator <- df_fin %>% 
#  dplyr::select(site, clim_class)

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

## Get summary tables ---------------------------------------------------------

# how many sites do not have any stem density present? not evemn mature treees?
# Count rows and calculate proportions
total_sites <- length(unique(df_fin$site))

df_fin %>%
  dplyr::filter(sum_stems_mature > 0) %>%
  dplyr::count(adv_delayed) %>%
  mutate(
    proportion = n / total_sites,
    proportion_percent = (n / total_sites) * 100
  )


# what are dominant species : from rIVI?
# Create a table of dominant_species with counts and proportions
species_dominance_rIVI <- df_fin %>%
  count(dominant_species) %>% # from rIVI
  mutate(
    #proportion = n / sum(n),
    proportion_percent = (n / sum(n)) * 100
  ) %>%
  arrange(desc(n))

# Display the result
species_dominance_rIVI


## Species composition: Tables --------------

df_stem_sp_sum <- stem_dens_species_long_cluster %>% 
  group_by(cluster, Species) %>% 
  summarise(sum_stem_density = sum(stem_density, na.rm = t)) %>% 
  dplyr::filter(sum_stem_density>0)



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


# Find the top species with more then 5% share 
top_species_overall <- species_composition_overall %>%
  arrange(desc(share)) %>%  
#  slice_head(n = 7) #%>%  # Select the top X species
  dplyr::filter(share > 5)# %>%  # select species with share > 5%

top_species_overall_vect <- top_species_overall$Species

(top_species_overall_vect)

# add to gpkg

# Create a new column in the sf object
df_fin_clim_clust_xy$dominant_species_grouped <- ifelse(
  df_fin_clim_clust_xy$dominant_species %in% top_species_overall_vect,
  df_fin_clim_clust_xy$dominant_species,  # Keep the species name if it's in the top species
  "other"  # Otherwise, assign 'other'
)
#st_write(df_fin_clim_clust_xy, "outData/xy_clim_cluster.gpkg", layer = "df_fin", driver = "GPKG", append = FALSE)



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


# Find the top species > 5% share 
top_species_layer <- species_composition_layer %>%
  group_by(VegType) %>% 
  arrange(desc(share)) %>%  
  #slice_head(n = 7) #%>%  # Select the top X species
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
# colors for map:
#'#B0B0B0' # forest?
#'#F7F7F7' # background?

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
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Black border around facets
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

my_species_levels <-  levels(df_stem_sp_sum_ordered$Species)

# Add a log-transformed column for sum_stem_density
df_stem_sp_sum_ordered <- df_stem_sp_sum_ordered %>%
  mutate(log_sum_stem_density = log10(sum_stem_density + 1))  # Adding 1 to avoid log(0)

# test with log values
p_stem_density_species <- df_stem_sp_sum_ordered %>%
  ggplot(aes(x = log_sum_stem_density, y = Species, group = Species)) +
  geom_density_ridges(aes(fill = Species), alpha = 0.8, color = 'NA') +
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
    size = 0.2,
    linewidth = 0.1,
    stroke = 0.2,
    position = position_nudge(y = 0.3)  # Adjust position slightly
  ) +
  theme_classic() +
  labs(title = "",
       x = "Stem density (log10)\n[n/ha]",
       y = "")  +
  scale_x_continuous(
    labels = math_format(10^.x)  # Format x-axis labels as 10^3, 10^4, etc.
  ) +
   theme(    axis.text = element_text(size = 8),    # Axis tick labels
             axis.title = element_text(size = 8),   # Axis titles
             
     #axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate the y-axis text (Species)
       # panel.border = element_rect(color = "black", fill = NA, size = 1),  # Black border around facets
        strip.background = element_rect(color = "black", linewidth = 1),  # Black border around facet labels
        panel.grid.major = element_line(color = "grey", linetype = "dotted"),  # Add grey dashed major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        legend.position = 'none'#,
     #   plot.margin = margin(t = 7, r = 5, b = 5, l = 3)
     )  # Adjust plot margins) #+  # Hide the legend


p_stem_density_species

ggsave(filename = 'outFigs/p_stem_density_ridge_log.png', 
       plot = p_stem_density_species, 
       width = 3, height = 3.5, dpi = 300, bg = 'white')


## Bar plot with error plot ---------------------------------------------------------


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
  mutate(species_VegType = paste(Species, VegType_acc, sep = '_')) %>% 
  mutate(VegType_acc = factor(VegType_acc, 
                              levels = c('mat', 'juv', 'sap'))) 





## Make a barplot for variation in stem density : try with shaded colors ---------------------------

# Make a color gradient: shaded within species
library(colorspace)   # For gradient and color manipulation

# Generate gradients for each species and VegType (sap, juv, mat)
species_shaded_colors <- unlist(lapply(species_colors, function(base_color) {
  gradient_n_pal(c(lighten(base_color, 0.6), base_color, darken(base_color, 0.4)))(c(0.33, 0.66, 1))
}))

# Create combined names for each Species-VegType combination
names(species_shaded_colors) <- c(
  "frex_sap", "frex_juv", "frex_mat",
  "piab_sap", "piab_juv", "piab_mat",
  "pisy_sap", "pisy_juv", "pisy_mat",
  "fasy_sap", "fasy_juv", "fasy_mat",
  "besp_sap", "besp_juv", "besp_mat",
  "acps_sap", "acps_juv", "acps_mat",
  "soau_sap", "soau_juv", "soau_mat",
  "quro_sap", "quro_juv", "quro_mat",
  "potr_sap", "potr_juv", "potr_mat",
  "abal_sap", "abal_juv", "abal_mat"
)


# frex_sap  frex_juv  frex_mat  piab_sap  piab_juv  piab_mat  pisy_sap  pisy_juv  pisy_mat  fasy_sap  fasy_juv  fasy_mat  besp_sap  besp_juv  besp_mat 
# "#C02C3E" "#9B0023" "#87001D" "#E85447" "#C8372D" "#A33833" "#F98260" "#E5643A" "#C65228" "#FFB97B" "#EEA150" "#CF852A" "#FFE49D" "#EED076" "#CCAE48" 
# acps_sap  acps_juv  acps_mat  soau_sap  soau_juv  soau_mat  quro_sap  quro_juv  quro_mat  potr_sap  potr_juv  potr_mat  abal_sap  abal_juv  abal_mat 
# "#DCF28F" "#C9DE79" "#A7BC53" "#ADE072" "#98CA5B" "#7CAC3C" "#72C86F" "#5FB05D" "#51944F" "#36A861" "#158E4A" "#0B793D" "#2A7D4D" "#016133" "#02532C" 



p_bar_IRQ_shaded <- stem_dens_species_long_cluster %>%  
  dplyr::filter(Species %in% top_species_overall_vect[1:7]) %>% 
  dplyr::filter(stem_density > 0) %>% 
  mutate(Species = factor(Species, levels = rev(my_species_levels ))) %>%  # Set ordertop_species_overall_vect[1:7]
  mutate(
   # Species = factor(Species, levels = top_species_overall_vect[1:7]),  # Set order
    species_veg = paste(Species, VegType_acc, sep = "_")  # Combine Species and VegType_acc
  ) %>%
  ggplot(aes(y = VegType_acc, x = stem_density/1000 , fill = species_veg)) +
  
  # Bar plot for median with error bars for IQR
  stat_summary(
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
    fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
    geom = "bar", 
    position = position_dodge2(width = 0.2, padding = 0.2),  # Dodge for horizontal bars
    color = "black",
    orientation = "y" ,
    linewidth = 0
    #width = 0.7
  ) +
  stat_summary(
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # Add error bars
    fun.max = function(x) quantile(x, 0.75), 
    geom = "pointrange", 
    position = position_dodge2(width = 0.2, padding = 0.2),
    orientation = "y",
  #   color = "black",  # Black outline for points
     shape = 21,  # Shape 21 is a circle with a fill and border
     size = 0.2,
    linewidth = 0.1,
  stroke = 0.2
  ) +
  # Adjust scales and themes
  facet_grid(Species ~ ., switch = "y") +  # Facet by species
  scale_fill_manual(values = species_shaded_colors ) +  # colors species_colors
  theme_classic() +
  theme(    axis.text = element_text(size = 6),    # Axis tick labels
            axis.title = element_text(size = 7),   # Axis titles
            
    legend.position = 'none',
    strip.background = element_blank(),  # Remove facet background
    strip.text.y.left = element_text(face = "plain", angle = 0, vjust = 1),  # Horizontal facet labels
    strip.placement = "outside",  # Place labels outside
    plot.margin = margin(t = 7, r = 5, b = 5, l = 3),  # Adjust plot margins
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  
  # Labels and formatting
  labs(x = "Stem density\n[*1000 n/ha]", y = "")

p_bar_IRQ_shaded

p_bar_IRQ_color_sp <- stem_dens_species_long_cluster %>%  
  dplyr::filter(Species %in% top_species_overall_vect[1:7]) %>% 
  dplyr::filter(stem_density > 0) %>% 
  mutate(Species = factor(Species, levels = my_species_levels)) %>%  # Set order
  mutate(
    Species = factor(Species, levels = top_species_overall_vect[1:7]),  # Set order
    species_veg = paste(Species, VegType_acc, sep = "_")  # Combine Species and VegType_acc
  ) %>%
  ggplot(aes(y = VegType_acc, x = stem_density , fill = Species)) +
  
  # Bar plot for median with error bars for IQR
  stat_summary(
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
    fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
    geom = "bar", 
    position = position_dodge2(padding = 0.2),  # Dodge for horizontal bars
    color = "black",
    orientation = "y" #,
    #width = 0.7
  ) +
  stat_summary(
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # Add error bars
    fun.max = function(x) quantile(x, 0.75), 
    geom = "pointrange", 
    position = position_dodge2(padding = 0.2),
    orientation = "y",
    color = "black",  # Black outline for points
    shape = 21,  # Shape 21 is a circle with a fill and border
    size = 0.2,
    linewidth = 0.1,
    stroke = 0.2
  ) +
  # Adjust scales and themes
  facet_grid(Species ~ ., switch = "y") +  # Facet by species
  scale_fill_manual(values = species_colors ) +  # colors species_colors
  theme_classic() +
  theme(
    axis.text = element_text(size = 8),    # Axis tick labels
    axis.title = element_text(size = 8),   # Axis titles
    
    legend.position = 'none',
    strip.background = element_blank(),  # Remove facet background
    strip.text.y.left = element_text(face = "bold", angle = 0, vjust = 1),  # Horizontal facet labels
    strip.placement = "outside",  # Place labels outside
    plot.margin = margin(t = 7, r = 5, b = 5, l = 7),  # Adjust plot margins
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  
  # Labels and formatting
  labs(x = "Stem density\n[n/ha]", y = "") +
  scale_x_continuous(labels = math_format(10^.x))  # Format x-axis labels





# export plot alternatives -------------------------------------------------


# export plot together with vertical layers
p_stem_dens_composition_old <- ggarrange(p_stem_density_species, 
                                         p_species_vert_layer,
                                         #align = 'hv', 
                                        # axis = "tb",
                                         labels = c("[a]", "[b]"), font.label = list(size = 8, face = "plain") )

p_stem_dens_composition_old

ggsave(filename = 'outFigs/p_stem_dens_composition.png', 
       plot = p_stem_dens_composition_old, 
       width = 5, height = 3.5, dpi = 300, bg = 'white')




# export plot together with vertical layers
p_stem_dens_composition_stem_dens_shaded <- ggarrange(p_stem_density_species, 
                                                      p_bar_IRQ_shaded ,
                                                     # align = 'hv', 
                                                     # axis = "lr",
                                                      labels = c("[a]", "[b]"), font.label = list(size = 8, face = "plain"),
                                                     widths = c(1, 1.5) )

p_stem_dens_composition_stem_dens_shaded

ggsave(filename = 'outFigs/p_stem_dens_composition_stem_dens_shaded.png', 
       plot = p_stem_dens_composition_stem_dens_shaded, 
       width = 4, height = 2.5, dpi = 300, bg = 'white')



# export plot together with vertical layers
p_stem_dens_composition_stem_dens_sp <- ggarrange(p_stem_density_species, 
                                                  p_bar_IRQ_color_sp ,
                                                  #align = 'hv', 
                                                  #axis = "tb",
                                                  labels = c("[a]", "[b]"), font.label = list(size = 8, face = "plain") )

p_stem_dens_composition_stem_dens_sp

ggsave(filename = 'outFigs/p_stem_dens_composition_stem_dens_sp.png', 
       plot = p_stem_dens_composition_stem_dens_sp, 
       width = 4.5, height = 3.5, dpi = 300, bg = 'white')







## Tableplot:  Stem density per species and vertical class: --------------------------------

# Get a summary table:
summary_stem_dens_VegType <- stem_dens_species_long_cluster %>%
  dplyr::filter(Species %in% top_species_overall_vect[1:7]) %>%
  dplyr::filter(stem_density > 0) %>%
  mutate(Species = factor(Species, levels = top_species_overall_vect[1:7])) %>%
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
summary_stem_dens_spec <- stem_dens_species_long_cluster %>%
  dplyr::filter(Species %in% top_species_overall_vect[1:7]) %>%
  dplyr::filter(stem_density > 0) %>%
  mutate(Species = factor(Species, levels = top_species_overall_vect[1:7])) %>%
  group_by(Species) %>%
  summarise(
    min = min(stem_density, na.rm =T),
    max = max(stem_density, na.rm =T),
    mean     = mean(stem_density, na.rm =T),
    sd = sd(stem_density, na.rm =T),
    Median = median(stem_density, na.rm =T),
    IQR = IQR(stem_density)#,
    
  ) %>%
  arrange(Species)

summary_stem_dens_spec


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
sjPlot::tab_df(summary_stem_dens_spec_formated, 
       file = "outTable/summary_stem_dens_spec.doc", 
       title = "Summary of Stem Density per Species",
       show.rownames = FALSE)







p_stem_density_error <- stem_dens_species_long_cluster %>%  
  dplyr::filter(Species %in% top_species_overall_vect[1:7]  ) %>% 
  dplyr::filter(stem_density > 0) %>% 
  mutate(Species = factor(Species, 
                          levels = top_species_overall_vect[1:7])) %>%
  
  ggplot(aes(x = stem_density, 
             y = VegType_acc, 
             fill = VegType_acc,
             color = VegType_acc,
             group = species_VegType )) +
  # geom_violin(trim = T) +
  # geom_density_ridges(aes(fill = species_VegType ), alpha = 0.5) +
  stat_summary(
    aes(x = stem_density), 
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
    fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
    geom = "pointrange", 
    #color = VegType, 
    size = 0.2#,
    #position = position_nudge(y = .2)  # Adjust position slightly
  ) +
  #scale_color_manual(values = colorRampPalette(brewer.pal(11, "RdYlGn"))(3)) +  # Apply the color palette based on seral type
  facet_grid(Species ~ ., switch = "y") +
  
  # Adjust theme
  theme_classic() +
  theme(
    legend.position = 'none',
    strip.background = element_blank(),  # Remove background from facet labels
    strip.text.y.left = element_text(face = "bold", angle = 0,vjust = 1),  # Make facet labels bold and horizontal
    strip.placement = "outside",  # Place facet labels further outside
    
    # Expand plot margins to allow space for labels on the left
    plot.margin = margin(t = 5, r = 5, b = 5, l = 7),  # Increase left margin
    
    # Ensure xy lines appear only on axes
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  labs(x = "", y = "")

p_stem_density_error
ggsave(filename = 'outFigs/p_stem_density_error.png', 
       plot = p_stem_density_error, 
       width = 3, height = 3.5, dpi = 300, bg = 'white')






## Density plot per vertical class: -----------------------------------------------

stem_dens_species_long_cluster %>%  
  dplyr::filter(Species %in% top_species_overall_vect[1:7]  ) %>% 
  dplyr::filter(stem_density > 0) %>% 
  mutate(Species = factor(Species, levels = top_species_overall_vect[1:7])) %>%  # Set order
  ggplot(aes(x = log_stem_density , y = VegType_acc, 
             fill = VegType_acc,
             #color = VegType_acc,
             group = species_VegType )) +
  geom_density_ridges(aes(fill = Species), alpha = 0.5, color = 'NA') +
  #scale_fill_manual(values = species_colors) +
  stat_summary(
    aes(x = log_stem_density), 
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
    fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
    geom = "pointrange", 
    size = 0.1,
    position = position_nudge(y = .2)  # Adjust position slightly
  ) +
  #scale_color_manual(values = colorRampPalette(brewer.pal(11, "RdYlGn"))(3)) +  # Apply the color palette based on seral type
  facet_grid(Species ~ ., switch = "y") +
  scale_x_continuous(
    labels = math_format(10^.x)  # Format x-axis labels as 10^3, 10^4, etc.
  ) +
  # Adjust theme
  theme_classic() +
  theme(
    legend.position = 'none',
    strip.background = element_blank(),  # Remove background from facet labels
    strip.text.y.left = element_text(face = "bold", angle = 0,vjust = 1),  # Make facet labels bold and horizontal
    strip.placement = "outside",  # Place facet labels further outside
    
    # Expand plot margins to allow space for labels on the left
    plot.margin = margin(t = 5, r = 5, b = 5, l = 7),  # Increase left margin
    
    # Ensure xy lines appear only on axes
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  labs(x = "", y = "")
 
  



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



#### Species composition : Richness ----------------------
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
 


#### Species composition: country -----------------------------------

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
  ungroup() %>% 
  dplyr::filter(share >5)


# Create a color palette with 12 colors based on RdYlGn
my_colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))(length(unique(species_composition$Species)))

# Create the stacked bar plot
p_species_distribution_country <- species_composition %>% 
  dplyr::filter(share >5) %>% 
  ggplot(
                                 aes(x = country, 
                                     y = share, 
                                     fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  geom_text(aes(label = ifelse(share >= 5, paste0(round(share, 1), "%"), "")),
            position = position_stack(vjust = 0.5),  # Labels inside the bars
            size = 3, color = "black") +  # Adjust text size and color
  labs(x = "", y = "Percentage", 
       fill = "Species",
       title = "") +
  scale_fill_manual(values = my_colors) +
  #scale_fill_manual(values = species_colors) +  # Apply the color palette based on seral type
  theme_classic() +  # Use a clean theme
  theme(
    legend.position = 'right',
    # axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) 

p_species_distribution_country


ggsave(filename = 'outFigs/fig_p_species_distribution_global_country.png', 
       plot = p_species_distribution_country, 
       width = 7, height = 5, dpi = 300, bg = 'white')




## Structure ---------------------------------------------------------

### Vertical structure  --------------------------------------

#### Get presence absence data for individual layers 

# Group by cluster and VegType, check if stem_density > 0
vert_class_presence_absence <- stem_dens_species_long_cluster %>%
  group_by(cluster, VegType_acc) %>%
  summarise(has_stem_density = sum(stem_density > 0, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  # Filter only where there is stem density present
  dplyr::filter(has_stem_density) %>%
  dplyr::select(cluster, VegType_acc) %>% 
  mutate(Presence = 1) %>%  # Add a column with 1 indicating presence
  pivot_wider(names_from = VegType_acc, values_from = Presence, values_fill = 0)  # Convert to wide format, filling NAs with 0



# Now find clusters where no VegType has stem_density > 0
# Find clusters where no vegType has stem_density
all_clusters <- stem_dens_species_long_cluster %>%
  ungroup() %>% 
  dplyr::select(cluster, VegType_acc) %>% 
  distinct(cluster)

# Combine clusters with and without stem density
vert_class_presence_absence_fin <- vert_class_presence_absence  %>% 
  right_join(all_clusters)



# visualize vertical classes by UpSEt plot

### Upset plot -----------------------------------

library(UpSetR)
library("ComplexUpset")


library(ggupset)

tidy_movies %>%
  distinct(title, year, length, .keep_all=TRUE) %>%
  #str()
  ggplot(aes(x=Genres)) +
  geom_bar() +
  scale_x_upset(n_intersections = 20)
  
vert_class_presence_absence_fin %>% 
  ggplot(aes(x=))

  
upset_data <- vert_class_presence_absence_fin %>%
  rowwise() %>%
  mutate(
    sets = list(
      c(
        if (sap == 1) "sap" else NULL,
        if (juv == 1) "juv" else NULL,
        if (mat == 1) "mat" else NULL
      )
    )
  ) %>%
  ungroup() %>%
  filter(lengths(sets) > 0)  # Remove rows with no sets


# full data:

# Step 2: Select only the columns with presence/absence data
upset_data <- vert_class_presence_absence_fin %>% dplyr::select(sap, juv, mat) %>% 
  as.data.frame()

# Step 3: Create the UpSet plot
p_upset_test <- upset(upset_data, sets = c('sap', 'juv', 'mat'), order.by = "freq")
p_upset_test

# Define your color palette with colorBrewer
upset_colors <- brewer.pal(3, "Set2")  # Change 'Set2' to your desired palette


# Create the UpSet plot with ComplexUpset
# Create the UpSet plot
p_upset_test <- ComplexUpset::upset(
  upset_data, 
  intersect = c('sap', 'juv', 'mat'),  # Specify your sets
  name = "Intersection", 
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(vjust = -0.5)  # Adjusts text position
    ) + scale_fill_manual(values = upset_colors)  # Use a manual fill scale based on color palette
  ),
  themes = ggplot2::theme_minimal() +  # Minimal theme for ggplot
    theme(legend.position = "none")     # Remove unnecessary legend if not needed
)

# Display the plot
p_upset_test

#### upset data test ---------------------------------------------------------

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








### Descriptive plots  -----------------------------------------------------------

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
                   prop = count/n_subplot *100) %>% 
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







# 5. Drivers ---------------------------------
# for all regeneration (juveniles and saplings)
# split table in two: drivers for advanced (> 1000 stems of juveniles/ha)
#                     drivers for delayed regeneration (<50 stems/ha: saplings  + juveniles)  

## Prep final table --------------

df_fin <- df_fin %>% 
  mutate(country_full    = factor(country_full),
         country_abbr    = factor(country_abbr),
         country_pooled  = factor(country_pooled),
         region          = factor(region))

# account for climatic variables - aggregated on 9 km resolution:
df_fin <- df_fin %>%
  mutate(clim_grid = factor(paste(round(tmp,3), round(prcp,3), sep = "_")))  # Create a unique identifier for temp/precip combinations



# Make sure df_fin has no missing values, if necessary
df_fin <- na.omit(df_fin)

# make  regeneration classes: advanced vs delayed
df_fin <- df_fin %>% 
  mutate(adv_delayed = factor(ifelse(stem_regeneration <= 50, "Delayed", 
                                     ifelse(sum_stems_juvenile >= 1000, "Advanced", "Other")),
                               levels = c("Delayed", "Other", "Advanced")))

# include mature trees: present/absent
df_fin$sum_stems_mature_pres_abs <- ifelse(df_fin$sum_stems_mature > 0, 1, 0)

# scaled mature trees 0-1
df_fin$sum_stems_mature_scaled <- (df_fin$sum_stems_mature - min(df_fin$sum_stems_mature)) / 
  (max(df_fin$sum_stems_mature) - min(df_fin$sum_stems_mature))


# check my categories?
df_fin %>% 
  dplyr::select(site, 
                stem_density,
                stem_regeneration,
                sum_stems_juvenile,
                sum_stems_sapling, 
                #sum_stems_mature,
                adv_delayed
               
  ) %>% 
  View()




#### check for multicollinearity -----------------------------------------------------

library(car)
selected_data <- df_fin %>%
  dplyr::select(stem_density, 
                spei12,
                tmp,  
                prcp 
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

# tw distribution is the best


# TW has a better fit, also can handle zero!
# select main variables as predictors 
predictor_vars_sub <- c(
  "spei12",
  "tmp", 
  "prcp", 
  
  # 
  "salvage_intensity",
  "protection_intensity",
  "management_intensity",
  "distance_edge", 
  
  # disturbance severity est
  "disturbance_severity", # from RS 
 
  # disturbance severity based on residual trees: 
   "mature_dist_severity",  # cover of residual trees over subplots
   "sum_stems_mature",       # stem density of mature trees
   "sum_stems_mature_scaled" , # mature trees stems: scaled
   #"sum_stems_mature_pres_abs",  # mature trees present/absent
  
  "clay_extract", 
  "depth_extract", 
   
   # site info
  "av.nitro",
  #"richness",
  #'rIVI',
  
    "n_vertical")



# run univariate models to find teh best predictors :

# Define response and predictor variables
response_var <- "stem_regeneration" # Replace with your actual response variable name
predictor_vars <- predictor_vars_sub

# Create an empty data frame to store results
results <- data.frame(
  Predictor = character(),
  AIC = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each predictor variable
for (predictor in predictor_vars) {
  # Create the formula
  formula <- as.formula(paste(response_var, "~ s(", predictor, ", k = 3)"))
  
  # Fit the GAM model with Tweedie distribution
  model <- gam(formula, family = tw(), data = df_fin)
  
  # Extract AIC
  aic <- AIC(model)
  
  # Store the results
  results <- rbind(results, data.frame(Predictor = predictor, AIC = aic))
}

# Sort results by AIC (lower is better)
results <- results[order(results$AIC), ]

# Display the results
print(results)



# correlation ;

# Initialize a data frame to store results
cor_results <- data.frame(
  Predictor = character(),
  Correlation = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each predictor variable
for (predictor in predictor_vars) {
  # Calculate correlation, handling NA values
  cor_value <- cor(df_fin[[response_var]], df_fin[[predictor]], use = "complete.obs", method = "spearman")
  
  # Store the results
  cor_results <- rbind(cor_results, data.frame(Predictor = predictor, Correlation = cor_value))
}

# Sort results by absolute correlation
cor_results <- cor_results[order(-abs(cor_results$Correlation)), ]

# Display the results
print(cor_results)

# check new predictors
# Subset the data with the variables of interest
eda_data <- df_fin[, c("stem_regeneration", "sum_stems_mature", "sum_stems_mature_scaled", "sum_stems_mature_pres_abs", "mature_dist_severity",
                       "disturbance_severity", "distance_edge")]

# Generate a scatterplot matrix
pairs(eda_data, main = "")



## Models: simplify table just for stem_regeneration  ---------------------------------------------------------------------------------
# test drivers: simplify the analysis:
# Subset the data

df_stem_regeneration2 <- df_fin %>% 
  dplyr::select(all_of(c("stem_regeneration", predictor_vars_sub,
                         "country_pooled", "clim_grid",  "x", "y")))


# Centering the variables in your data frame
# Centering all relevant continuous variables in your data frame
df_stem_regeneration2 <- df_stem_regeneration2 %>%
  mutate(
    prcp_c = prcp - mean(prcp, na.rm = TRUE),
    tmp_c = tmp - mean(tmp, na.rm = TRUE),
    spei12_c = spei12 - mean(spei12, na.rm = TRUE),
    distance_edge_c = distance_edge - mean(distance_edge, na.rm = TRUE),
    disturbance_severity_c = disturbance_severity - mean(disturbance_severity, na.rm = TRUE),
    mature_dist_severity_c = mature_dist_severity - mean(mature_dist_severity, na.rm = TRUE),
    clay_extract_c = clay_extract - mean(clay_extract, na.rm = TRUE),
   # av.nitro_c = av.nitro - mean(av.nitro, na.rm = TRUE),
    depth_extract_c = depth_extract - mean(depth_extract, na.rm = TRUE)#,
    #management_intensity_c = management_intensity - mean(management_intensity, na.rm = TRUE)
  )

summary(df_stem_regeneration2)

## Drivers: regeneration density pooled ---------------------------------------------


# Check the structure of the data
str(df_stem_regeneration2)

# Check for missing values
sum(is.na(df_stem_regeneration2))

# Inspect the distribution of each predictor
hist(df_stem_regeneration2$prcp, main="Histogram of Precipitation (prcp)", xlab="prcp")
hist(df_stem_regeneration2$tmp, main="Histogram of Temperature (tmp)", xlab="tmp")

# check correlation between predictors
library(corrplot)

# Select the relevant predictors from your data frame
predictors <- df_stem_regeneration2 %>%
  dplyr::select(prcp, tmp, spei12, distance_edge, depth_extract, 
                disturbance_severity, mature_dist_severity , 
                sum_stems_mature ,
                depth_extract , clay_extract, av.nitro, management_intensity)

# Calculate the correlation matrix
correlation_matrix <- cor(predictors, use = "complete.obs")

# Display the correlation matrix
print(correlation_matrix)

# check correlations between disturbance severty from EU map and site-base (from presence/absence of mature trees)
cor(df_fin$disturbance_severity, df_fin$sum_stems_mature, method = "spearman")
cor(df_fin$disturbance_severity, df_fin$sum_stems_mature, method = "kendall") # better if many values lies in same tieghts
cor(df_fin$disturbance_severity, df_fin$mature_dist_severity, method = "kendall") # better if many values lies in same tieghts
plot(df_fin$disturbance_severity, df_fin$mature_dist_severity)

par(mfrow = c(1, 2))
plot(df_fin$prcp, df_fin$clim_grid)
plot(df_fin$tmp, df_fin$clim_grid )
dev.off()

ggplot(df_stem_regeneration2, aes(x = clay_extract, y = stem_regeneration     )) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue")


# get stats: distnace to edge and severity typesper country
ggplot(df_stem_regeneration2, aes(x = country_pooled , y = distance_edge     )) +
  geom_boxplot()
  
# get stats: distnace to edge and severity typesper country
ggplot(df_stem_regeneration2, aes(x = country_pooled , y = disturbance_severity     )) +
  geom_boxplot()

# get stats: distnace to edge and severity typesper country
ggplot(df_stem_regeneration2, aes(x = country_pooled , y = mature_dist_severity     )) +
  geom_boxplot()



# split analysis in two ??? nope! keep random effects to have both scales in ----------------------------------------------------------
# one for climate: on grid 9x9, for disturbance severity vs distance_edge - on local basis

df_median <- df_stem_regeneration2 %>%
  group_by(clim_grid) %>%
  summarise(
    # Compute the median for numeric columns
    across(where(is.numeric), \(x) median(x, na.rm = TRUE)),
    # Retain the first value for factor columns
    across(where(is.factor), ~ first(.)),
    # Retain the first value for character columns
    across(where(is.character), ~ first(.))
  )

# View the result
head(df_median)

# averaged values lead to completely different results! importance of tmp, not prec
m_med_int <- gam(stem_regeneration ~
               s(prcp, k = 5) + s(tmp, k = 5) + 
               s(spei12, k = 5) + 
               #s(distance_edge, k = 5) +
               s(depth_extract, k = 4) +
               #s(disturbance_severity, k =5) +
               #s(mature_dist_severity, k = 5) +
               s(clay_extract, k = 5) +
               s(av.nitro, k =5) +
               ti(tmp, prcp, k = 5) +
               #ti(disturbance_severity, distance_edge, k = 10) +
               #ti(disturbance_severity_c, prcp_c, k = 5) +
               s(management_intensity,by = country_pooled, k = 4) + 
               s(country_pooled, bs = "re") + 
               s(x,y) #+
               #s(clim_grid, bs = "re") 
             ,
             family = tw(), method = "REML", data = df_median)


plot(m_med_int, page = 1)


# test between large scale an local drivers:  keep them in separate model

# averaged values lead to completely different results! importance of tmp, not prec
m_macro_int <- gam(stem_regeneration ~
                   s(prcp, k = 5) + s(tmp, k = 5) + 
                   s(spei12, k = 5) + 
                   #s(distance_edge, k = 5) +
                   s(depth_extract, k = 4) +
                   #s(disturbance_severity, k =5) +
                   #s(mature_dist_severity, k = 5) +
                   s(clay_extract, k = 5) +
                   s(av.nitro, k =5) +
                   ti(tmp, prcp, k = 5) +
                   #ti(disturbance_severity, distance_edge, k = 10) +
                   #ti(disturbance_severity_c, prcp_c, k = 5) +
                   s(management_intensity,by = country_pooled, k = 4) + 
                   s(country_pooled, bs = "re") + 
                   s(x,y) #+
                 #s(clim_grid, bs = "re") 
                 ,
                 family = tw(), method = "REML", data = df_stem_regeneration2)


# # averaged values lead to completely different results! importance of tmp, not prec
m_local_int <- gam(stem_regeneration ~
                     #s(prcp, k = 5) + s(tmp, k = 5) + 
                     #s(spei12, k = 5) + 
                     s(distance_edge, k = 5) +
                     #s(depth_extract, k = 4) +
                     s(disturbance_severity, k =5) +
                     s(mature_dist_severity, k = 5) +
                     #s(clay_extract, k = 5) +
                     #s(av.nitro, k =5) +
                     #ti(tmp, prcp, k = 5) +
                     ti(disturbance_severity, distance_edge, k = 10) +
                     #ti(disturbance_severity_c, prcp_c, k = 5) +
                     s(management_intensity,by = country_pooled, k = 4) + 
                     s(country_pooled, bs = "re") + 
                     s(x,y) #+
                   #s(clim_grid, bs = "re") 
                   ,
                   family = tw(), method = "REML", data = df_stem_regeneration2)


m_local_int2 <- gam(stem_regeneration ~
                     #s(prcp, k = 5) + s(tmp, k = 5) + 
                     #s(spei12, k = 5) + 
                     s(distance_edge, k = 5) +
                     #s(depth_extract, k = 4) +
                     s(disturbance_severity, k =5) +
                     s(mature_dist_severity, k = 5) +
                     #s(clay_extract, k = 5) +
                     #s(av.nitro, k =5) +
                     #ti(tmp, prcp, k = 5) +
                     ti(mature_dist_severity, distance_edge, k = 5) +
                     #ti(disturbance_severity_c, prcp_c, k = 5) +
                     s(management_intensity,by = country_pooled, k = 4) + 
                     s(country_pooled, bs = "re") + 
                     s(x,y) #+
                   #s(clim_grid, bs = "re") 
                   ,
                   family = tw(), method = "REML", data = df_stem_regeneration2)



# test across scales:

macro_model <- gam(
  stem_regeneration ~ s(prcp, k = 5) + s(tmp, k = 5) + s(clim_grid, bs = "re"),
  family = tw(),
  data = df_stem_regeneration2
)

micro_model <- gam(
  stem_regeneration ~ s(disturbance_severity, k = 5) + s(distance_edge, k = 5),
  family = tw(),
  data = df_stem_regeneration2
)

AIC(macro_model, micro_model, m_combined)
summary(macro_model)
summary(micro_model)


m_combined <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(disturbance_severity, k = 5) + s(distance_edge, k = 5) +
    ti(prcp, disturbance_severity, k = 5) +  # Interaction term
    ti(distance_edge, disturbance_severity, k = 5) +  # Interaction term
    ti(prcp, tmp, k = 5) +  # Interaction term
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  data = df_stem_regeneration2
)

# test with diferent types of disturbance severity -------------------
m_combined_mature <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(mature_dist_severity , k = 5) + s(distance_edge, k = 5) +
    ti(prcp, mature_dist_severity , k = 5) +  # Interaction term across cales
    ti(distance_edge, mature_dist_severity , k = 5) +  # Interaction term
    ti(prcp, tmp, k = 5) +  # Interaction term
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  data = df_stem_regeneration2
)

m_combined_mature_basic <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(mature_dist_severity , k = 5) + s(distance_edge, k = 5) +
   # ti(prcp, mature_dist_severity , k = 5) +  # Interaction term
    ti(distance_edge, mature_dist_severity , k = 5) +  # Interaction term
    ti(prcp, tmp, k = 5) +  # Interaction term
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  data = df_stem_regeneration2
)

# !!!
m_combined_mature_manag <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(mature_dist_severity , k = 5) + s(distance_edge, k = 5) +
    ti(prcp, mature_dist_severity , k = 5) +  # Interaction term across cales
    ti(distance_edge, mature_dist_severity , k = 5) +  # Interaction term
    ti(prcp, tmp, k = 5) +  # Interaction term
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  data = df_stem_regeneration2
)


m_combined_mature_manag2 <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(mature_dist_severity , k = 5) + s(distance_edge, k = 5) +
    ti(prcp, mature_dist_severity , k = 5) +  # Interaction term across cales
    ti(distance_edge, mature_dist_severity , k = 3) +  # Interaction term
    ti(prcp, tmp, k = 10) +  # Interaction term
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  data = df_stem_regeneration2
)

m_combined_mature_dist <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(mature_dist_severity , k = 5) + s(distance_edge, k = 5) +
    ti(prcp, mature_dist_severity , k = 5) +  # Interaction term across cales
    ti(prcp, distance_edge , k = 5) +  # Interaction term across cales
    ti(distance_edge, mature_dist_severity , k = 5) +  # Interaction term
    ti(prcp, tmp, k = 5) +  # Interaction term
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  data = df_stem_regeneration2
)

m_combined_mature_stems <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(sum_stems_mature  , k = 5) + s(distance_edge, k = 5) +
    ti(prcp, sum_stems_mature  , k = 5) +  # Interaction term
    ti(distance_edge, sum_stems_mature  , k = 5) +  # Interaction term
    ti(prcp, tmp, k = 5) +  # Interaction term
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  data = df_stem_regeneration2
) 
# !!!

m_combined_mature_stems_sc <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(sum_stems_mature_scaled  , k = 5) + s(distance_edge, k = 5) +
    ti(prcp, sum_stems_mature_scaled  , k = 5) +  # Interaction term
    ti(distance_edge, sum_stems_mature_scaled  , k = 5) +  # Interaction term
    ti(prcp, tmp, k = 5) +  # Interaction term
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  data = df_stem_regeneration2
)


AIC(m_combined_mature_manag,
    m_combined_mature_stems, m_combined_mature, m_combined,m_combined_mature_stems_sc ,m_combined_mature_basic, m_combined_mature_dist)




table(df_fin$adv_delayed)

m <- m_combined_mature_manag
summary(m)
vis.gam(m, view = c("prcp", "tmp"))

vis.gam(m, view = c("prcp", "mature_dist_severity"))

vis.gam(m, view = c("distance_edge", "mature_dist_severity"))
plot(m, page = 1)

df_pred_test <- ggpredict(m, #int_re_fin_dist, 
                          terms = c(#"tmp_c", 'prcp_c' 
                            #'prcp',
                            #'tmp[8,9,10]'
                            #'mature_dist_severity',
                            #'distance_edge[100,200,300]'#,
                            #'disturbance_severity'
                            'distance_edge',
                            "mature_dist_severity[0.4,0.7,0.9]"
                            #'mature_dist_severity[0.2,0.5,0.7]'#,
                            
                            
                            #tmp,
                            #"distance_edge",
                            #"disturbance_severity[0.7,0.8,0.9]",
                            #"disturbance_severity", 
                            #        "prcp"
                          ))

ggplot(df_pred_test, aes(x = x, y = predicted )) +
  geom_line(aes(color = group), linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) #+
  #facet_grid(.~group    )

ggplot(df_pred_test, aes(x = x, y = predicted )) +
  geom_line(aes(color = group), linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2)# +
facet_grid(.~group    )



# store the best model for regeneration density
fin.m.reg.density <- int_re_fin 

vis.gam(fin.m.reg.density, view = c("prcp", "tmp"), plot.type = "persp",
        main = "Interaction between Precipitation and Temperature",
        zlab = "Stem Regeneration", xlab = "Precipitation", ylab = "Temperature")


# chack variability withing clim grid

# Calculate standard deviation of stem_regeneration for each clim_grid
stem_variability <- df_stem_regeneration2 %>%
  group_by(clim_grid) %>%
  summarise(
    sd_stem_regeneration = sd(stem_regeneration, na.rm = TRUE)
  )

# View variability
head(stem_variability)

# Plot variability
ggplot(stem_variability, aes(x = clim_grid, y = sd_stem_regeneration)) +
  geom_bar(stat = "identity") +
  labs(title = "Variability in Stem Regeneration by Climate Grid",
       x = "Climate Grid",
       y = "SD of Stem Regeneration") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# quick test 11/12/2024 -------------
# remove disturbence perc interaction, run with raw values, not a _c
m_simple <- gam(stem_regeneration ~
                    s(prcp, k = 5) + s(tmp, k = 5) + 
                     s(spei12, k = 5) + 
                    s(distance_edge, k = 5) +
                    s(depth_extract, k = 4) +
                    s(disturbance_severity, k =5) +
                  s(mature_dist_severity, k = 5) +
                  s(sum_stems_mature, k = 5) +
                    s(clay_extract, k = 5) +
                    s(av.nitro, k =5) +
                    #ti(tmp, prcp, k = 5) +
                    #ti(disturbance_severity_c, distance_edge_c, k = 10) +
                    #ti(disturbance_severity_c, prcp_c, k = 5) +
                    #s(management_intensity,by = country_pooled, k = 4) + 
                    s(country_pooled, bs = "re")# + 
                    #s(x,y) +
                    #s(clim_grid, bs = "re") 
                  ,
                  family = tw(), method = "REML", data = df_stem_regeneration2)

appraise(m_simple)
summary(m_simple)
k.check(m_simple)

m_int <- gam(stem_regeneration ~
               s(prcp, k = 5) + s(tmp, k = 5) + 
               s(spei12, k = 5) + 
               s(distance_edge, k = 5) +
               s(depth_extract, k = 4) +
               s(disturbance_severity, k =5) +
               s(mature_dist_severity, k = 5) +
               s(clay_extract, k = 5) +
               s(av.nitro, k =5) +
               ti(tmp, prcp, k = 5) +
               ti(disturbance_severity, distance_edge, k = 10) +
               #ti(disturbance_severity_c, prcp_c, k = 5) +
               s(management_intensity,by = country_pooled, k = 4) + 
               s(country_pooled, bs = "re") + 
             s(x,y) +
             s(clim_grid, bs = "re") 
             ,
             family = tw(), method = "REML", data = df_stem_regeneration2)


m_int_mature <- gam(stem_regeneration ~
               s(prcp, k = 5) + s(tmp, k = 5) + 
               s(spei12, k = 5) + 
               s(distance_edge, k = 5) +
               s(depth_extract, k = 4) +
               s(disturbance_severity, k =5) +
               s(mature_dist_severity, k = 5) +
               s(clay_extract, k = 5) +
               s(av.nitro, k =5) +
               ti(tmp, prcp, k = 5) +
               ti(mature_dist_severity, distance_edge, k = 5) +
               #ti(disturbance_severity_c, prcp_c, k = 5) +
               s(management_intensity,by = country_pooled, k = 4) + 
               s(country_pooled, bs = "re") + 
               s(x,y) +
               s(clim_grid, bs = "re") 
             ,
             family = tw(), method = "REML", data = df_stem_regeneration2)
AIC(m_int_mature, m_int)

summary(m_int)


# check the role of clim_grid
m_clim_grid <- gam(stem_regeneration ~
                      # s(prcp, k = 5) + s(tmp, k = 5) + 
                      # s(spei12, k = 5) + 
                      # s(distance_edge, k = 5) +
                      # s(depth_extract, k = 4) +
                      # s(disturbance_severity, k =5) +
                      # s(mature_dist_severity, k = 5) +
                      # s(clay_extract, k = 5) +
                      # s(av.nitro, k =5) +
                      # ti(tmp, prcp, k = 5) +
                      # ti(mature_dist_severity, distance_edge, k = 5) +
                      # #ti(disturbance_severity_c, prcp_c, k = 5) +
                      # s(management_intensity,by = country_pooled, k = 4) + 
                      # s(country_pooled, bs = "re") + 
                      # s(x,y) +
                      s(clim_grid, bs = "re") 
                    ,
                    family = tw(), method = "REML", data = df_stem_regeneration2)




m <- m_clim_grid
plot(m, page = 1)
vis.gam(m)
k.check(m)
gam.check(m)
summary(m)


# old ---------------


# remove disturbence perc interaction, run with raw values, not a _c
int_re_fin <- gam(stem_regeneration ~
                                 s(prcp, k = 15) + s(tmp, k = 15) + 
                                # s(spei12_c, k = 4) + 
                                 s(distance_edge, k = 15) +
                                 #s(depth_extract_c, k = 4) +
                                 s(disturbance_severity, k = 15) +
                                 s(clay_extract, k = 5) +
                                 #s(av.nitro_c, k =5) +
                                 ti(tmp, prcp, k = 5) +
                                 #ti(disturbance_severity_c, distance_edge_c, k = 10) +
                                 #ti(disturbance_severity_c, prcp_c, k = 5) +
                                 s(management_intensity,by = country_pooled, k = 4) + 
                                 s(country_pooled, bs = "re") + 
                                 s(x,y) +
                                 s(clim_grid, bs = "re") 
                               ,
                               family = tw(), method = "REML", data = df_stem_regeneration2)


int_re_fin_dist <- gam(stem_regeneration ~
                    s(prcp, k = 15) + s(tmp, k = 15) + 
                     s(spei12, k = 4) + 
                    s(distance_edge, k = 15) +
                    #s(depth_extract_c, k = 4) +
                    s(disturbance_severity, k = 15) +
                    s(clay_extract, k = 5) +
                    #s(av.nitro_c, k =5) +
                    ti(tmp, prcp, k = 5) +
                    ti(disturbance_severity, distance_edge, k = 10) +
                    #ti(disturbance_severity_c, prcp_c, k = 5) +
                    s(management_intensity,by = country_pooled, k = 4) + 
                    s(country_pooled, bs = "re") + 
                    s(x,y) +
                    s(clim_grid, bs = "re") 
                  ,
                  family = tw(), method = "REML", data = df_stem_regeneration2)

AIC(int_re_fin_dist, int_re_fin_mature, int_re_fin_mature_c, int_re_fin_mature_low5)
# calculate severity from mature trees:
int_re_fin_mature <- gam(stem_regeneration ~
                         s(prcp, k = 15) + s(tmp, k = 15) + 
                         s(spei12, k = 4) + 
                         s(distance_edge, k = 15) +
                         #s(depth_extract_c, k = 4) +
                         s(mature_dist_severity, k = 5) +
                        s(clay_extract, k = 5) +
                         #s(av.nitro_c, k =5) +
                         ti(tmp, prcp, k = 5) +
                         ti(mature_dist_severity, distance_edge, k = 5) +
                         #ti(disturbance_severity_c, prcp_c, k = 5) +
                         s(management_intensity,by = country_pooled, k = 5) + 
                         s(country_pooled, bs = "re") + 
                         s(x,y) +
                         s(clim_grid, bs = "re") 
                       ,
                       family = tw(), method = "REML", data = df_stem_regeneration2)


# calculate severity from mature trees: - low k
int_re_fin_mature_low5 <- gam(stem_regeneration ~
                           s(prcp, k = 5) + s(tmp, k = 5) + 
                           s(spei12, k = 4) + 
                           s(distance_edge, k = 5) +
                           #s(depth_extract_c, k = 4) +
                           s(mature_dist_severity, k = 5) +
                           s(clay_extract, k = 5) +
                           #s(av.nitro_c, k =5) +
                           ti(tmp, prcp, k = 5) +
                           ti(mature_dist_severity, distance_edge, k = 5) +
                           #ti(disturbance_severity_c, prcp_c, k = 5) +
                           s(management_intensity,by = country_pooled, k = 5) + 
                           s(country_pooled, bs = "re") + 
                           s(x,y) +
                           s(clim_grid, bs = "re") 
                         ,
                         family = tw(), method = "REML", data = df_stem_regeneration2)


summary(int_re_fin_mature_low5)
# use centered values:
# calculate severity from mature trees:
int_re_fin_mature_c <- gam(stem_regeneration ~
                           s(prcp_c, k = 15) + s(tmp_c, k = 15) + 
                           # s(spei12_c, k = 4) + 
                           s(distance_edge, k = 15) +
                           #s(depth_extract_c, k = 4) +
                           s(mature_disturbance_severity_c, k = 5) +
                           # s(clay_extract, k = 5) +
                           #s(av.nitro_c, k =5) +
                           ti(tmp_c, prcp_c, k = 5) +
                           ti(mature_disturbance_severity_c, distance_edge_c, k = 5) +
                           #ti(disturbance_severity_c, prcp_c, k = 5) +
                           s(management_intensity_c,by = country_pooled, k = 5) + 
                           s(country_pooled, bs = "re") + 
                           s(x,y) +
                           s(clim_grid, bs = "re") 
                         ,
                         family = tw(), method = "REML", data = df_stem_regeneration2)

AIC(int_re_fin_mature_c, int_re_fin_mature,int_re_fin_dist)

# 



int_re_fin_spei <- gam(stem_regeneration ~
                    s(prcp, k = 15) + s(tmp, k = 15) + 
                     s(spei12_c, k = 4) + 
                    s(distance_edge, k = 15) +
                    #s(depth_extract_c, k = 4) +
                    s(disturbance_severity, k = 15) +
                    s(clay_extract, k = 5) +
                    #s(av.nitro_c, k =5) +
                    ti(tmp, prcp, k = 5) +
                    #ti(disturbance_severity_c, distance_edge_c, k = 10) +
                    #ti(disturbance_severity_c, prcp_c, k = 5) +
                    s(management_intensity,by = country_pooled, k = 4) + 
                    s(country_pooled, bs = "re") + 
                    s(x,y) +
                    s(clim_grid, bs = "re") 
                  ,
                  family = tw(), method = "REML", data = df_stem_regeneration2)

AIC(int_re_fin_spei,int_re_fin)



AIC(int_re_fin, int_re_dist_prcp_simpl4 , int_re_dist_prcp_simpl3, int_re_dist_prcp_simpl2, int_re_dist_prcp_simpl)
k.check(int_re_dist_prcp_simpl4)
summary(int_re_dist_prcp_simpl4)

# quick plotting --------------------------------------

table(df_fin$adv_delayed)

m <- int_re_fin_mature_low5
vis.gam(m)
plot(m, page = 1)

df_pred_test <- ggpredict(m, #int_re_fin_dist, 
                          terms = c(#"tmp_c", 'prcp_c' 
                            'distance_edge',
                            'mature_dist_severity[0.4,0.6]'#,
                           
                                    
                            #tmp,
                            #"distance_edge",
                            #"disturbance_severity[0.7,0.8,0.9]",
                            #"disturbance_severity", 
                            #        "prcp"
                            ))

ggplot(df_pred_test, aes(x = x, y = predicted )) +
  geom_line(aes(color = group), linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2)# +
  facet_grid(.~group    )

  ggplot(df_pred_test, aes(x = x, y = predicted )) +
    geom_line(aes(color = group), linewidth = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2)# +
  facet_grid(.~group    )
  
  

# store the best model for regeneration density
fin.m.reg.density <- int_re_fin 

vis.gam(fin.m.reg.density, view = c("prcp", "tmp"), plot.type = "persp",
        main = "Interaction between Precipitation and Temperature",
        zlab = "Stem Regeneration", xlab = "Precipitation", ylab = "Temperature")

# # ti(disturbance_severity, distance_edge): Interaction between disturbance severity and distance to edge
# vis.gam(fin.m.reg.density, view = c("disturbance_severity_c", "distance_edge_c"), plot.type = "persp",
#         zlab = "Stem Regeneration")

# ti(disturbance_severity, distance_edge): Interaction between disturbance severity and distance to edge
vis.gam(fin.m.reg.density, view = c("disturbance_severity_c", "prcp_c"), plot.type = "persp",
        zlab = "Stem Regeneration")


# check again correlation between predictors:

# Load necessary libraries
library(corrplot)

# Select predictors from your data frame
predictors <- df_stem_regeneration2[, c("prcp", "tmp",# "spei12_c", 
                                        "distance_edge", "depth_extract",
                                        "disturbance_severity", 
                                        "clay_extract", "av.nitro", 
                                        "management_intensity")]

# Calculate correlation matrix
correlation_matrix <- cor(predictors,  method = "spearman", use = "complete.obs")
windows()
# Plot the correlation matrix
corrplot(correlation_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, 
         title = "Correlation Matrix of Predictors",
         addCoef.col = "black", number.cex = 0.7)

# Export as a PNG
png("outFigs/correlation_matrix_plot.png", width = 800, height = 800, res = 300)
plot.new()
corrplot(correlation_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, 
         title = "Correlation Matrix of Predictors",
         addCoef.col = "black", number.cex = 0.7)
dev.off()#summary(interaction_model_4)


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


# Plot  ---------------------------------------------------------------------------
y_lab = 'Stem density [#/ha]'

summary(fin.m.reg.density)

# show only 95% quatile for stem density
# Define the quantiles for stem_density and tmp columns
quantiles_stem_density <- quantile(df_stem_regeneration2$stem_regeneration  , 
                                   probs = c(0, 0.95), na.rm = TRUE)
#quantiles_tmp <- quantile(df_stem_regeneration2$tmp, probs = c(0, 0.95), na.rm = TRUE)

# Filter the DataFrame to keep rows within these quantile ranges
filtered_df_plot <- df_stem_regeneration2 %>%
  dplyr::filter(stem_regeneration     >= quantiles_stem_density[1] & stem_regeneration  <= quantiles_stem_density[2])
# ,tmp >= quantiles_tmp[1] & tmp <= quantiles_tmp[2]

# Display the filtered data
filtered_df

# Create plots for individual variables
plot_prcp <- create_plot(fin.m.reg.density, "prcp", 
                         df_stem_regeneration2,
                         #filtered_df_plot, 
                         "precipitation 0.005 **", line_color = "blue", fill_color = "blue")
plot_prcp
plot_tmp <- create_plot(fin.m.reg.density, "tmp", df_stem_regeneration2, 
                        "temperature 0.002 **", line_color = "red", fill_color = "red")
plot_disturbance_severity <- create_plot(fin.m.reg.density, 
                                         "disturbance_severity", df_stem_regeneration2, 
                                         "dist severity 0.02*", line_color = "purple", fill_color = "purple")
plot_distance_edge <- create_plot(fin.m.reg.density, "distance_edge", 
                                  df_stem_regeneration2, "dist edge  0.17 ns",  line_color = "darkgreen", fill_color = "darkgreen")


# Create interaction plots
plot_interaction1 <- create_interaction_plot(fin.m.reg.density, 
                                             c("prcp", "tmp[8,9,10]"), "prcp & tmp 0.04 *", 
                                             df_stem_regeneration2) +
  labs(color = "tmp", fill = "tmp") 
(plot_interaction1)




library(cowplot)
combined_plot <- plot_grid(plot_prcp, plot_tmp,plot_interaction1,
                        #    plot_spei12, 
                        plot_disturbance_severity, plot_distance_edge, 
                        #    plot_management_intensity, 
                        ncol = 2,nrow = 2,
                        labels = c("[a]","[b]", "[c]","[d]", "[e]"), 
                        label_size = 8, label_fontface = "plain")

combined_plot


# Save the combined plot
ggsave('outFigs/fig_regen_pool_drivers.png', plot = combined_plot, 
       width = 6, height = 6.5, bg = 'white')



# Identify random effects using the model's "smooth" component
smooth_terms <- summary(fin.m.reg.density)$s.table

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
sjPlot::tab_model(fin.m.reg.density,
                  show.re.var = TRUE,        # Show the variance components
                  #show.icc = TRUE,           # Show Intraclass Correlation Coefficient
                  #show.dev = TRUE,           # Show deviance
                  pred.labels = c("Intercept", pred_labels), # Replace smooth term labels
                  dv.labels = paste0("Explained Deviance: ", round(100 * summary(fin.m.reg.density)$dev.expl, 2), "%"), 
                  file = "outTable/full_drivers_reg_stem_density.doc")











## Wilcox: Boxplot sites differences: delayed vs advaced ----------------------------

# two categories: count how many plots i have?

#prop.table(table(df_fin$adv_delayed_wider))
prop.table(table(df_fin$adv_delayed))
table(df_fin$adv_delayed)


# Select relevant columns including 'RegenerationStatus' and the desired variables
variables_to_plot <- c(
  "tmp",
  "prcp",     
  "distance_edge" , 
  "disturbance_severity" 
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
    IQR = IQR(Value, na.rm = TRUE)) %>% 
  arrange(adv_delayed)
print(summary_stats_narrow)

  df_fin %>%
  na.omit() %>% 
  dplyr::select(all_of(c(variables_to_plot, 'adv_delayed'))) %>% 
  gather(key = "Variable", value = "Value", -adv_delayed) %>%
  group_by(Variable) %>%
  summarise(
    Median = median(Value, na.rm = TRUE),
    IQR = IQR(Value, na.rm = TRUE)) 

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

# get indicators for delayed/advanced: ------------------------------------------
df_delayed_advanced <- df_fin %>% 
  dplyr::select(site, country, adv_delayed)

fwrite(df_delayed_advanced, 'outTable/df_delayed_advanced.csv')

# summarize stem deregeneration, rIVI, richness, n_vert
df_fin %>% 
  group_by(adv_delayed) %>% 
  summarise(med_stem_density = median(stem_density),
            iqr_stem_density = IQR(stem_density),
            med_rIVI = median(rIVI),
            iqr_rIVI = IQR(rIVI),
            med_richness = median(richness),
            iqr_richness = IQR(richness),
            med_n_vertical = median(n_vertical),
            iqr_n_vertical = IQR(n_vertical)
            )


### test differences between sites: --------------------------------------
# Step 1: Reshape the data to long format
df_long_narrow <- df_fin %>%
  na.omit() %>% 
  dplyr::select(tmp, 
                prcp, 
                #spei12,
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
  labs(title = "",
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


# filter out outliers
df_long_narrow_filtered <- df_long_narrow %>%
  group_by(adv_delayed, Variable) %>%
  dplyr::filter(
    Value > quantile(Value, 0.25) - 1.5 * IQR(Value) & 
      Value < quantile(Value, 0.75) + 1.5 * IQR(Value)
  ) %>%
  ungroup()


df_fin %>% 
  dplyr::select(site, sum_stems_juvenile, sum_stems_sapling, 
                sum_stems_mature, stem_regeneration, adv_delayed) %>% 
  dplyr::filter(adv_delayed == "Delayed") %>% 
  View()
  #count(adv_delayed)

p_boxplot_wilcox_outliers <- 
 # df_long_narrow_filtered %>% 
  df_long_narrow %>% 
  mutate(Variable = factor(Variable, 
                           levels = c('prcp', 'tmp', "distance_edge", "disturbance_severity"))) %>% 
  ggboxplot(
          x = "adv_delayed", y = "Value", 
          fill = "adv_delayed", 
          palette = c("#A50026", 
                      "#FDAE61",
                      "#006837"),
          alpha = 0.5,
          facet.by = "Variable", 
          scales = "free_y", 
          ylab = "Values", xlab = "Regeneration Status",
          outlier.shape = NA,  # Hide outliers
          size = 0.2) +
    geom_jitter(size = 0.1, alpha = 0.5, aes(color = adv_delayed)) +
  scale_color_manual(values = c("#A50026", 
                                "#FDAE61",
                                "#006837")) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", 
                     size = 2,
                     label.x = 1.5) +  # Position labels between the groups
  labs(title = "",
       x = "", y = "Vals") +
  theme(
    legend.position = 'none',
    text = element_text(size = 5),         # Set all text to size 3
    axis.text = element_text(size = 5),    # Axis tick labels
    axis.title = element_text(size = 5),   # Axis titles
    strip.text = element_text(size = 5),   # Facet labels
    legend.text = element_text(size = 5),  # Legend text
    plot.title = element_text(size = 5),   # Plot title
    strip.background = element_blank(),    # Remove the box around facet names
    strip.placement = "outside",           # Optional: Move facet label outside the plot area
    #panel.border = element_blank(),        # Remove top and right border
    panel.grid = element_blank(),          # Optional: Remove grid lines for a cleaner look
    axis.line = element_line(color = "black", linewidth = 0.5)  # Add bottom and left axis lines
  )
 #"#A50026" "#DA362A" "#F46D43" "#FDAE61" "#FEE08B" "#D9EF8B" "#A6D96A" "#66BD63" 
# Save the combined plot (optional)



# filtered outliers
p_boxplot_wilcox_rm_outliers <- 
   df_long_narrow_filtered %>% 
  #df_long_narrow %>% 
  mutate(Variable = factor(Variable, 
                           levels = c('prcp', 'tmp', "distance_edge", "disturbance_severity"))) %>% 
  ggboxplot(
    x = "adv_delayed", y = "Value", 
    fill = "adv_delayed", 
    palette = c("#A50026", 
                "#FDAE61",
                "#006837"),
    alpha = 0.5,
    facet.by = "Variable", 
    scales = "free_y", 
    ylab = "Values", xlab = "Regeneration Status",
    outlier.shape = NA,  # Hide outliers
    size = 0.2) +
  geom_jitter(size = 0.1, alpha = 0.5, aes(color = adv_delayed)) +
  scale_color_manual(values = c("#A50026", 
                                "#FDAE61",
                                "#006837")) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", 
                     size = 2,
                     label.x = 1.5) +  # Position labels between the groups
  labs(title = "",
       x = "", y = "Vals") +
  theme(
    legend.position = 'none',
    text = element_text(size = 5),         # Set all text to size 3
    axis.text = element_text(size = 5),    # Axis tick labels
    axis.title = element_text(size = 5),   # Axis titles
    strip.text = element_text(size = 5),   # Facet labels
    legend.text = element_text(size = 5),  # Legend text
    plot.title = element_text(size = 5),   # Plot title
    strip.background = element_blank(),    # Remove the box around facet names
    strip.placement = "outside",           # Optional: Move facet label outside the plot area
    #panel.border = element_blank(),        # Remove top and right border
    panel.grid = element_blank(),          # Optional: Remove grid lines for a cleaner look
    axis.line = element_line(color = "black", linewidth = 0.5)  # Add bottom and left axis lines
  )
#"#A50026" "#DA362A" "#F46D43" "#FDAE61" "#FEE08B" "#D9EF8B" "#A6D96A" "#66BD63" 
# Save the combined plot (optional)
# Save the plot ensuring text sizes are preserved


# Save the plot ensuring text sizes are preserved
ggsave("outFigs/p_boxplot_wilcox_with_outliers.png", plot = p_boxplot_wilcox_outliers , 
       width = 3, height = 3.2, units = "in", dpi = 300, 
       bg = 'white', scale = 1)



# Save the plot ensuring text sizes are preserved
ggsave("outFigs/p_boxplot_wilcox_no_outliers.png", plot = p_boxplot_wilcox_rm_outliers, 
       width = 3, height = 3.2, units = "in", dpi = 300, 
       bg = 'white', scale = 1)


#### Wilcox tablle -----------------------------------------------

# Perform Wilcoxon pairwise comparisons and extract p-values
wilcox_results <- df_long_narrow %>%
  #df_long_narrow_filtered %>% 
  group_by(Variable) %>%
  summarise(
    p_values = list(pairwise.wilcox.test(Value, adv_delayed, p.adjust.method = "bonferroni")$p.value)
  )

# Convert the Wilcoxon results into a tidy format
tidy_wilcox_results <- wilcox_results %>%
  rowwise() %>%
  mutate(
    tidy_p_values = list(as.data.frame(as.table(p_values)))
  ) %>%
  unnest(tidy_p_values) %>%
  dplyr::rename(Group1 = Var1, Group2 = Var2, p_value = Freq) %>%
  arrange(Variable, Group1, Group2)

# Display the tidy results
tidy_wilcox_results


df_fin %>% 
  group_by(adv_delayed) %>% 
  summarize(mean = mean(stem_density),
            sd = sd(stem_density),
            median = median(stem_density),
            IQR = IQR(stem_density))


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
save(fin.m.reg.density, df_stem_regeneration2,
     file = "outData/stem_density_models.RData")

# 4. future developmenet -------------------------------------------

# compare current species composition with future ones form wessely

# 
head(stem_dens_species_long_cluster)

# get species merging table
species_look_up_simple<- read.csv("rawData/tree_sp_simple.csv", sep = ';')

# get species merging table
species_look_up_full<- read.csv("rawData/tree_sp_field_wessely_merged.csv", sep = ';')


# identify what species are present per plot
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


# get me an overview: what are the species that are mostly lost?
# Find species not suitable under any RCP but present currently
# Count how often each species is out of climatic suitability
species_out_of_suitability_count <- df_compare_future_species %>%
  ungroup() %>%
  dplyr::filter(current == 1 & rcp26 == 0 & rcp45 == 0 & rcp85 == 0) %>%
  dplyr::count(acc, name = "count") %>%
  mutate(proportion = count / total_sites*100) %>%
  arrange(desc(count))

# Display the results
species_out_of_suitability_count


# which species will remain?
# Calculate the presence proportion for each scenario
species_presence_proportion <- df_compare_future_species %>%
  ungroup() %>%
  dplyr::filter(current == 1) %>%
  dplyr::group_by(acc) %>%
  dplyr::summarise(
    current_count = sum(current == 1, na.rm = T),
    rcp26_count = sum(rcp26 == 1, na.rm = T),
    rcp45_count = sum(rcp45 == 1, na.rm = T),
    rcp85_count = sum(rcp85 == 1, na.rm = T),
    current_proportion = (current_count / total_sites) * 100,
    rcp26_proportion = (rcp26_count / total_sites) * 100,
    rcp45_proportion = (rcp45_count / total_sites) * 100,
    rcp85_proportion = (rcp85_count / total_sites) * 100
  ) %>%
  arrange(desc(current_count))

# Display the results
View(species_presence_proportion)


# Identify plots with no species presence under each RCP scenario
#plots_no_species <- 
  df_compare_future_species %>%
  group_by(site) %>%
  summarise(
    current_species_present = sum(current == 1, na.rm = TRUE),
    rcp26_species_present = sum(rcp26 == 1, na.rm = TRUE),
    rcp45_species_present = sum(rcp45 == 1, na.rm = TRUE),
    rcp85_species_present = sum(rcp85 == 1, na.rm = TRUE)
  ) %>%
  summarise(
    no_species_rcp26 = sum(rcp26_species_present == 0),
    no_species_rcp45 = sum(rcp45_species_present == 0),
    no_species_rcp85 = sum(rcp85_species_present == 0)
  )

# Display the number of plots with no species under each scenario
plots_no_species



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
