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
library(cowplot)

# Make a color gradient: shaded within species
library(colorspace)   # For gradient and color manipulation




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


# update labels
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


# Print the color assignments for confirmation
print(species_colors)

top_species_site_share <- c("piab", "fasy", "quro", "pisy", "soau", "acps", "potr", "abal", "besp", "lade")

# Assign colors to each species in the order of `top_species_site_share$Species`
species_colors <- setNames(
  reversed_colors,
  c("piab", "fasy", "quro", "pisy", "soau", "acps", "potr", "abal", "besp", "lade")
)



# Read data -----------------------------------------------------------------------

# get vegetation data
load("outData/veg.Rdata")

# final tables on site level
df_fin <- fread('outData/indicators_for_cluster_analysis.csv')

# add manually created regions
regions_manual <- terra::vect('rawData/regions_manual2.gpkg')
regions_manual_df <- as.data.frame(regions_manual)

df_fin <- df_fin %>% 
  full_join(regions_manual_df) %>% 
  rename(region_manual = name)

# ad disturbance severity based on presence/absence of mature trees 
#df_mature_dist_severity <- fread('outData/disturb_severity_mature.csv')

# Create the opposite vector: listing residual trees
#df_mature_dist_severity$residual_mature_trees <- 1 - df_mature_dist_severity$mature_dist_severity


#df_fin <- df_fin %>% 
#  left_join(df_mature_dist_severity, by = c('site' = 'cluster')) %>% 
#  mutate(time_since_disturbance = 2023 - disturbance_year)

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
#print(df_fin_clim_clust_xy)

# Step 3: Export the data as a GeoPackage (GPKG)
#st_write(df_fin_clim_clust_xy, "outData/xy_clim_cluster.gpkg", layer = "df_fin", driver = "GPKG", append = FALSE)

# get only a ataframe of teh climate clusters and sites for easy merging to detiailed veg data
#clim_cluster_indicator <- df_fin %>% 
#  dplyr::select(site, clim_class)

#fwrite(df_fin, 'outTable/df_fin.csv')

prop.table(table(df_fin$disturbance_year))





## add to gpkg: 10 prevailing species 
# 
# # Create a new column in the sf object
# df_fin_clim_clust_xy$dominant_species_grouped <- ifelse(
#   df_fin_clim_clust_xy$dominant_species %in% top_species_site_share$Species,
#   df_fin_clim_clust_xy$dominant_species,  # Keep the species name if it's in the top species
#   "other"  # Otherwise, assign 'other'
# )
#st_write(df_fin_clim_clust_xy, "outData/xy_clim_cluster.gpkg", layer = "df_fin", driver = "GPKG", append = FALSE)


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


### rIVI -----------------------
#### rIVI How many plots per dominant tree species ? (rIVI) -----------
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


#### get plots share by most frequent dominant species: rIVI -------------------------------
# Filter species with > 5% proportion

# Create the barplot
p_bar_species_dominance <- species_dominance_rIVI %>%
  dplyr::filter(proportion_percent > 3) %>% 
  mutate(dominant_species = reorder(dominant_species, -proportion_percent)) %>%  # Reorder by descending proportion
  ggplot(aes(x = 1, y = proportion_percent, fill = dominant_species)) +
  geom_bar(stat = "identity", color = 'black') + 
  scale_fill_manual(values = species_colors) +  # Apply custom color palette
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  geom_text(aes(label = paste0(round(proportion_percent, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 2) +  # Add labels with share percentages
  labs(
    title = "",
    fill = 'Dominant tree\nspecies (rIVI)',
    x = NULL,
    y = "Share of plots [%]"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis text since it's a single bar
        axis.title.x = element_blank(),
        
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8, face = "plain")
  )


## Species richness --------------------------------------------

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
  summarise(site_count = n_distinct(site)) %>%
  mutate(proportion = site_count / sum(site_count) * 100)  # Calculate proportions

richness_summary

# Create a stacked bar plot
p_bar_richness_groups <- ggplot(richness_summary, 
                                aes(x = 1, 
                                    y = proportion, 
                                    fill = richness_category)) +
  geom_bar(stat = "identity", color ="black" 
             ) +
  geom_text(aes(label = paste0(round(proportion, 1), "%", '(', site_count, ')')), 
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
  summarise(site_count = n_distinct(site)) %>%
  mutate(proportion = site_count / sum(site_count) * 100)  # Calculate proportions

# Create a stacked bar plot
p_bar_vertical_groups <- ggplot(vert_str_summary, aes(x = 1, y = proportion, 
                                                      fill = vert_category)) +
  geom_bar(stat = "identity", color = "black"
             ) +
  geom_text(aes(label = paste0(round(proportion, 1), "%", '(', site_count, ')')), 
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


ggarrange(p_bar_vertical_groups,p_bar_species_dominance, p_bar_richness_groups,  
          labels = c("[a]", "[b]", "[c]"),  # Add labels
          ncol = 3,
          font.label = list(size = 10, face = "plain"))

## Species composition: Tables --------------

df_stem_sp_sum <- stem_dens_species_long_cluster %>% 
  group_by(cluster, Species) %>% 
  # summarize across vertical groups
  summarise(sum_stem_density = sum(stem_density, na.rm = t)) #%>% 
  #dplyr::filter(sum_stem_density>0)

# do i still have sites that have no stem density? 

df_stem_sp_sum %>% 
  dplyr::filter(sum_stem_density==0)

### Share of plots by species -----------------------------------------------------

# Count the unique clusters (sites) for each species
species_site_share <- df_stem_sp_sum %>%
  dplyr::filter(sum_stem_density >0) %>% 

  group_by(Species) %>%
  summarise(site_count = n_distinct(cluster)) %>%
  ungroup()


# Add a column for the share (percentage) of sites
species_site_share <- species_site_share %>%
   mutate(share_of_sites = (site_count / total_sites) * 100)

# Sort the results in descending order of share
species_site_share <- species_site_share %>%
  arrange(desc(share_of_sites))

# Print the result
print(species_site_share)

top_species_site_share <- species_site_share %>% 
  arrange(desc(share_of_sites)) %>%
  slice_head(n = 10) # Select the top 7 rows



## Create a horizontal bar plot -------------------------------------------------
p_sites_share_species <- species_site_share %>%
  arrange(desc(share_of_sites)) %>%
  slice_head(n = 10) %>% # Select the top X rows
  mutate(Species = factor(Species, levels = unique(Species))) %>% 
  ggplot(aes(x = share_of_sites, y = reorder(Species, share_of_sites))) +
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

p_sites_share_species



## Density plot of stem density  ----------------------------------

# Calculate median for each species and reorder the factor levels
df_stem_sp_sum_ordered <- df_stem_sp_sum %>%
  dplyr::filter(sum_stem_density >0) %>% 
  dplyr::filter(Species %in% top_species_site_share$Species ) %>%  #top_species_overall
  dplyr::group_by(Species) %>%
  dplyr::mutate(median_stem_density = median(sum_stem_density, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>%
  mutate(Species = factor(Species, levels = rev(top_species_site_share$Species))) # Set custom order
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
    x = expression("\nStem density(log"[10]*") [n ha"^-1*"]"),
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

# Print the plot
print(p_stem_density_species)


# % of recorded stems of juveniles and saplings ---------------------
stem_dens_species_long_cluster %>% 
  group_by(VegType) %>% 
  summarise(sum_vegType = sum(stem_density, na.rm = T)) %>% 
  mutate(total_stem_density = sum(sum_vegType),
         share = sum_vegType/total_stem_density)
# 
# # A tibble: 3 × 4
# VegType   sum_vegType total_stem_density  share
# <fct>           <dbl>              <dbl>  <dbl>
#   1 Mature         101000            6172500 0.0164
# 2 Juveniles      859125            6172500 0.139 
# 3 Saplings      5212375            6172500 0.844 

# get quantiles for stem density: here i need to account for empty clusters!
median(df_fin$stem_density)
mean(df_fin$stem_density)
#quantile(df_fin$stem_density)
quantile(df_fin$stem_density, probs = c(0, 0.05,0.1, 0.25, 0.5, 0.75,0.9,  0.95, 1), na.rm = TRUE)


## get % of stems per saplings/juveniles -----------------------------------------
# Summarize the total stem density per species 
species_composition_sapl_juv <- stem_dens_species_long_cluster %>%
  dplyr::filter(Species %in% top_species_site_share$Species) %>% 
  dplyr::filter(VegType != "Mature") %>% 
  group_by(cluster, Species,VegType) %>%
  summarize(sum_regeneration = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the stem density er regeneration class and share of each species
species_composition_sapl_juv <- 
  species_composition_sapl_juv %>%
  group_by(cluster) %>%
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
                          levels = rev(top_species_site_share$Species))) # Set custom order

  


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
    Species = factor(Species, levels = rev(top_species_site_share$Species)) # Custom species order
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

# View the results
print(overall_med_share)



## adjust colors -----------------------------------------------------------------

# Generate separate colors for Juveniles (darker) and Saplings (lighter)
#species_colors_juveniles <- species_colors

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
           label = "Sapling", color = "gray30", size = 3, fontface = "plain") + # Add Saplings label
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
p1 <- p_sites_share_species +
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

# Print combined plot
print(combined_plot)

# Save the combined plot as an image
ggsave(
  filename = "outFigs/combined_plot2.png",    # File name (change extension for different formats, e.g., .pdf)
  plot = combined_plot,             # Plot object to save
  width = 7,                       # Width of the saved plot in inches
  height = 3.5,                       # Height of the saved plot in inches
  dpi = 300,                        # Resolution (dots per inch)
  units = "in"                      # Units for width and height
)

ggsave(
  filename = "outFigs/combined_plot2.svg",  # Save as SVG
  plot = combined_plot,
  width = 7,
  height = 3.5,
  dpi = 300,
  units = "in"
)

# to properly align x axis labels: need to open in svg, and align them then esport to pdf and convert online to png! otherwise, 
# vertical lines between plots appeares (which i can't get rid of) 


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


# example for stem density calculation: 

dd <- data.frame(cluster = rep(1:3, each = 3),
                 species = rep(LETTERS[1:3], 3),
                 stem_density = c(3,2,0,
                                  0,0,0,
                                  10,0,0))

# include zeros or not?? -------------------------------------------------------------

# Summary Table
# Metric                                 | Include Zeros? | When to Use

# Total stem density per cluster         | Yes (✅)        | Landscape-level biomass/productivity
# Mean density per species (all clusters)| Yes (✅)        | How common/abundant is species overall?
# Mean density per species (when present)| No  (❌)        | How dense is species when it's present?
# Species richness (per cluster)         | No  (❌)        | Count species with stem_density > 0





## Make a barplot for variation in stem density (if species is present): per vertical classes  ---------------------------
p_bar_IRQ_color_sp <- stem_dens_species_long_cluster %>%  
  dplyr::filter(Species %in% top_species_site_share$Species) %>% 
  dplyr::filter(stem_density > 0) %>% 
  mutate(Species = factor(Species, levels = my_species_levels)) %>%  # Set order
  mutate(
   Species = factor(Species, levels = top_species_site_share$Species),  # Set order
   species_veg = paste(Species, VegType_acc, sep = "_") , # Combine Species and VegType_acc
   Species_label = factor(Species, levels = names(species_labels), labels = species_labels )
   ) %>%
 # mutate(Species = factor(Species, levels = names(species_labels), labels = species_labels)) %>% 

  ggplot(aes(y = VegType_acc, x = stem_density , fill = Species)) +
  
  # Bar plot for median with error bars for IQR
  stat_summary(
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
    fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
    geom = "bar", 
    position = position_dodge2(padding = 0.2),  # Dodge for horizontal bars
    color = "black"       ,
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
    color = "black"       ,  # Black outline for points
    shape = 21,  # Shape 21 is a circle with a fill and border
    size = 0.2,
    linewidth = 0.1,
    stroke = 0.2
  ) +
  # Adjust scales and themes
  facet_grid(Species_label ~ ., switch = "y") +  # Facet by species
  scale_fill_manual(values = species_colors ) +  # colors species_colors
  theme_classic() +
  theme(
    axis.text = element_text(size = 8),    # Axis tick labels
    axis.title = element_text(size = 8),   # Axis titles
    
    legend.position = 'none',
    strip.background = element_blank(),  # Remove facet background
    strip.text.y.left = element_text(face = "italic", 
                                     angle = 0, vjust = 1.3),  # Horizontal facet labels
    strip.placement = "outside",  # Place labels outside
    plot.margin = margin(t = 7, r = 5, b = 5, l = 7),  # Adjust plot margins
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) +
  # Labels and formatting
  labs(x = "Stem density\n[n/ha]", y = "")# +
 # scale_x_continuous(labels = math_format(10^.x))  # Format x-axis labels


p_bar_IRQ_color_sp


ggsave(
  filename = "outFigs/top10_sp_dens_vert_class.png",  
  plot = p_bar_IRQ_color_sp,
  width = 3.5,
  height = 3.5,
  dpi = 300,
  units = "in"
)


## Tableplot:  Stem density per species and vertical class (): --------------------------------

# Get a summary table if species is present:
summary_stem_dens_VegType <- stem_dens_species_long_cluster %>%
  dplyr::filter(Species %in% top_species_site_share$Species) %>%
  dplyr::filter(stem_density > 0) %>%
  mutate(Species = factor(Species, levels = top_species_site_share$Species)) %>%
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
  dplyr::filter(Species %in% top_species_site_share$Species) %>%
  dplyr::filter(stem_density > 0) %>%
  mutate(Species = factor(Species, levels = top_species_site_share$Species)) %>%
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

summary_stem_dens_spec %>% 
  arrange()


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
  mutate(total_stems = sum(sum_stems),  # Total stem density in each climate class
         share = (sum_stems / total_stems) * 100) %>%  # Calculate percentage share
  ungroup()




#### Species composition: country -----------------------------------
# get inication per inividual country 
df_fin_country_ind <- df_fin %>% 
  dplyr::select(site, country_abbr)

# Summarize the total stem density per species for each climate class
species_composition_country <- 
  stem_dens_species_long_cluster %>%
  dplyr::rename(site = cluster) %>% 
  #  head()
  left_join(df_fin_country_ind, by = c("site")) %>%
  #  View()
  group_by(Species, country_abbr) %>%
  summarize(sum_stems = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the total stem density per climate class and the share of each species
species_composition_country <- species_composition_country %>%
  group_by(country_abbr) %>%
  mutate(total_stems = sum(sum_stems),  # Total stem density in each climate class
         share = (sum_stems / total_stems) * 100) %>%  # Calculate percentage share
  ungroup() %>% 
  dplyr::filter(share >5) %>%
  group_by(country_abbr) %>% 
  arrange(share) %>% 
  mutate(Species = factor(Species))


# Create a color palette with 12 colors based on RdYlGn
my_colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))(length(unique(species_composition$Species)))
length(my_colors)

# Define species full Latin names
species_labels_full <- c(
  "abal" = "Abies alba",
  "acps" = "Acer pseudoplatanus",
  "besp" = "Betula sp.",
  "fasy" = "Fagus sylvatica",
  "piab" = "Picea abies",
  "pisy" = "Pinus sylvestris",
  "potr" = "Populus tremula",
  "prav" = "Prunus avium",
  "soau" = "Sorbus aucuparia",
  "frex" = "Fraxinus excelsior",
  "quro" = "Quercus robur/petraea",
  "lade" = "Larix decidua",
  "acpl" = "Acer platanoides",
  "cabe" = "Carpinus betulus"
)


# Create the stacked bar plot
p_species_distribution_country <- species_composition_country %>% 
  dplyr::filter(share >5) %>% 
  ggplot(    aes(x = country_abbr, 
                                     y = share, 
                                     fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  geom_text(aes(label = ifelse(share >= 5, paste0(round(share, 1), ""), "")),
            position = position_stack(vjust = 0.5),  # Labels inside the bars
            size = 3, 
            color ="black" 
              ) +  # Adjust text size and color
  labs(x = "", y = "Share [%]", 
       fill = "",
       title = "") +
  scale_fill_manual(values = rev(my_colors),
                    labels = species_labels_full) +  # Use full Latin names) +
  #scale_fill_manual(values = species_colors) +  # Apply the color palette based on seral type
  theme_classic() +  # Use a clean theme
  theme(
    legend.position = "bottom",  # Move legend to the bottom
    legend.text = element_text(size = 8, face = "italic"),  # Make legend text italic
    
    # axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  guides(fill = guide_legend(ncol = 3))  # Split legend into 3 columns

p_species_distribution_country


ggsave(filename = 'outFigs/fig_p_species_distribution_country_full.png', 
       plot = p_species_distribution_country, 
       width = 7, height = 6, dpi = 300, bg = 'white')




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
  

# full data:

# Step 2: Select only the columns with presence/absence data
upset_data <- vert_class_presence_absence_fin %>% 
  dplyr::select(sap, juv, mat) %>% 
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




# 5. Drivers ---------------------------------
# for all regeneration (juveniles and saplings)
# split table in two: drivers for advanced (> 1000 stems of juveniles/ha)
#                     drivers for delayed regeneration (<50 stems/ha: saplings  + juveniles)  

## Prep final table --------------

df_fin <- df_fin %>% 
  mutate(country_full    = factor(country_full),
         country_abbr    = factor(country_abbr),
         country_pooled  = factor(country_pooled),
         region          = factor(region),
         region_manual   = factor(region_manual)
         )

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
#df_fin$sum_stems_mature_pres_abs <- ifelse(df_fin$sum_stems_mature > 0, 1, 0)




# TW has a better fit, also can handle zero!
# select main variables as predictors 
predictor_vars_sub <- c(
  
  "spei1",
  "spei3",
  "spei12",
  "tmp", 
  "prcp", 
  
  "drought_tmp",
  "drought_prcp",
  
  "tmp_z", 
  "prcp_z", 
  
  # get drorught indiced : spei in 2018-2019
  "drought_spei1",
  "drought_spei3",
  "drought_spei12",
  
  # 
 # "salvage_intensity",
 # "protection_intensity",
#  "management_intensity",

  "distance_edge", 
  # disturbance severity est
  "disturbance_severity", # from RS 
 
  # disturbance severity based on residual trees: 
   #"mature_dist_severity",  # cover of residual trees over subplots
   #"sum_stems_mature",       # stem density of mature trees
   #"sum_stems_mature_scaled" , # mature trees stems: scaled
  # "residual_mature_trees",        
   #"sum_stems_mature_pres_abs",  # mature trees present/absent
  
  # site info
  "sand_extract",
  "clay_extract", 
  "depth_extract", 
  "av.nitro",
  
 # "time_since_disturbance",
  
  # anomalies per growing season for tmp and prcp
  "mean_grw_anm_prcp" ,        
  "mean_grw_anm_tmp"  ,        
  "max_grw_anm_prcp" ,        
  "max_grw_anm_tmp",
  "sd_grw_anm_prcp",           
  "sd_grw_anm_tmp" ,           


  # seasonality: CV - over year
  "cv_t2m",
  "cv_tp" #,
   
  #"richness",
  #'rIVI',
  
  #  "n_vertical"
  )






### run univariate models to find teh best predictors : -----------------------------

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
  model <- gam(formula, family = tw(), data = df_fin)
  
  # Extract AIC
  aic <- AIC(model)
  
  # Store the results
  AIC_results_univariate <- rbind(AIC_results_univariate, data.frame(Predictor = predictor, AIC = aic))
}

# Sort results by AIC (lower is better)
AIC_results_univariate <- AIC_results_univariate[order(AIC_results_univariate$AIC), ]

# Display the results
View(AIC_results_univariate)

plot(df_fin$drought_spei1  , df_fin$prcp)
cor(df_fin$drought_spei1  , df_fin$prcp, method = 'spearman')  # corr was small (0.12) between the spei12 and prcp
# high correlation (0.7) for drought_spei1& drougt_spei12 and prcp
 
sjPlot::tab_df(AIC_results_univariate,
               show.rownames = FALSE,
               file="outTable/univariate_AIC.doc",
               digits = 1) 



library(corrplot)
cor_matrix <- cor(df_fin[, ..predictor_vars], use = "complete.obs")
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7)



## Models: simplify table just for stem_regeneration  ---------------------------------------------------------------------------------
# test drivers: simplify the analysis:
# Subset the data

df_stem_regeneration2 <- df_fin %>% 
  dplyr::select(all_of(c("site", "stem_regeneration", predictor_vars_sub,
                         
                         "country_pooled","region", "region_manual", "clim_grid",  "x", "y")))


summary(df_stem_regeneration2)


#### check for multicollinearity -----------------------------------------------------

library(car)


df_small <- df_stem_regeneration2 %>% 
  dplyr::select(stem_regeneration, tmp, spei1, sd_grw_anm_tmp, #drought_prcp, 
                tmp_z, prcp_z)

lm_small <- lm(stem_regeneration ~ ., data = df_small)
vif(lm_small)



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
 # dplyr::select(where(is.numeric))
  dplyr::select(prcp, 
                tmp, 
                #spei12,
                drought_prcp, 
                drought_tmp, 
                #drought_spei12, 
                #drought_spei1,
                drought_spei3,
                sd_grw_anm_prcp,
                sd_grw_anm_tmp,
                #tmp_z, 
                #prcp_z,
                #cv_t2m,
                #cv_tp,
                #distance_edge, 
                depth_extract,
                disturbance_severity, #mature_dist_severity ,
                #sum_stems_mature ,
                depth_extract , 
                clay_extract, 
                av.nitro, #management_intensity
                )

# Calculate the correlation matrix
correlation_matrix <- cor(predictors, use = "complete.obs")

# Plot using corrplot
corrplot(correlation_matrix, method = "color", type = "upper",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         #addCoef.col = 
         , 
         tl.col = 
           , tl.srt = 45)

par(mfrow = c(1, 2))
plot(df_fin$prcp, df_fin$clim_grid)
plot(df_fin$tmp, df_fin$clim_grid )
dev.off()



# Export as a PNG
png("outFigs/correlation_matrix_plot.png", width = 800, height = 800, res = 300)
plot.new()
corrplot(correlation_matrix, method = "color", type = "upper",
         tl.col = 
           , tl.srt = 45, 
         title = "Correlation Matrix of Predictors",
         addCoef.col = 
           , number.cex = 0.7)
dev.off()#summary(interaction_model_4)






### keep only teh most meaningful model: 03192025 -----------
m_int_sev_edge_full_te_comb <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

# add soil predictors
m_int_sev_edge_full_te_comb_soil <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2

)

m_int_sev_edge_full_te_comb_soil1 <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    #te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m_soil_protect <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    #te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(protection_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m_soil_protect_int <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(protection_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)





m_fixed_soil <-  gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    #te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(protection_intensity,by = country_pooled, k = 4) + 
    #s(country_pooled, bs = "re") +
    #s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


AIC(m_int_sev_edge_full_te_comb_soil,m_int_sev_edge_full_te_comb_soil1, 
    m_int_sev_edge_full_te_comb, m_fixed_soil, 
    m_soil_protect,
    m_soil_protect_int)


AIC(m_fixed_soil,m_soil_protect )
summary(m_fixed_soil)
summary(m_soil_protect)



# 20250401 - completely remove management, also protection and salvage intensity ------------

m_rnd <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    #te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    #s(protection_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
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
    #s(protection_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
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
    #s(protection_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
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
    #s(protection_intensity,by = country_pooled, k = 4) + 
    #s(country_pooled, bs = "re") +
    #s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


# complete altarnative based on univariate AIC ---------
m_alt <- gam(
  stem_regeneration ~ 
    s(drought_prcp, k = 5) + s(tmp_z, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(drought_prcp,tmp_z, k = 5 ) +
    #s(protection_intensity,by = country_pooled, k = 4) + 
    #s(country_pooled, bs = "re") +
    #s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


AIC(m_rnd_ti_fixed, m_rnd,m_rnd_te, m_rnd_ti, m_alt)
summary(m_rnd)
summary(m_rnd_te)
summary(m_rnd_ti)
summary(m_rnd_ti_fixed)

# 20250603 add tmp anomalies into model: ----------------
# add tmp_z
m_rnd_ti_fixed_anm <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(tmp_z, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    #s(protection_intensity,by = country_pooled, k = 4) + 
    #s(country_pooled, bs = "re") +
    #s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

# add sd_deviation of teh growth
m_rnd_ti_fixed_anm2 <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(sd_grw_anm_tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    #s(protection_intensity,by = country_pooled, k = 4) + 
    #s(country_pooled, bs = "re") +
    #s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

# add drought_prcp
m_rnd_ti_fixed_anm3 <- gam(
  stem_regeneration ~ 
    s(drought_prcp, k = 5) + s(tmp, k = 5) +
    s(sd_grw_anm_tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    #s(protection_intensity,by = country_pooled, k = 4) + 
    #s(country_pooled, bs = "re") +
    #s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

# add drought_prcp
m_rnd_ti_fixed_anm4 <- gam(
  stem_regeneration ~ 
   # s(prcp, k = 5) +
    s(drought_prcp, k = 5) + s(tmp, k = 5) +
    s(sd_grw_anm_tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    #s(protection_intensity,by = country_pooled, k = 4) + 
    #s(country_pooled, bs = "re") +
    #s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


m_rnd_ti_fixed_anm5 <- gam(
  stem_regeneration ~ 
    # s(prcp, k = 5) +
    s(drought_prcp, k = 5) + s(tmp, k = 5) +
    s(sd_grw_anm_tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(drought_prcp,sd_grw_anm_tmp, k = 5 ) +
    #s(protection_intensity,by = country_pooled, k = 4) + 
    #s(country_pooled, bs = "re") +
    #s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

# drought_prcp

AIC( m_rnd_ti_fixed_anm, m_rnd_ti_fixed, m_rnd_ti_fixed_anm2,m_rnd_ti_fixed_anm3, m_rnd_ti_fixed_anm4, m_rnd_ti_fixed_anm5)
summary(m_rnd_ti_fixed_anm5)

plot(m_rnd_ti_fixed_anm5, pages = 1, scheme = 1, se = TRUE, rug = TRUE)
draw(m_rnd_ti_fixed_anm5)

# Basic 3D surface plot
vis.gam(m_rnd_ti_fixed_anm5,
        view = c("drought_prcp", "sd_grw_anm_tmp"),
        plot.type = "persp",  # use "contour" for flat plot
        color = "terrain",
        main = "Interaction: Drought × Temp Instability",
        zlab = "Effect on stem density",
        ticktype = "detailed")

vis.gam(m_rnd_ti_fixed_anm5,
        view = c("drought_prcp", "sd_grw_anm_tmp"),
        plot.type = "contour",
        color = "topo",
        main = "Interaction: Drought × Temp Instability")


# save teh best model
fin.m.reg.density <- m_rnd_ti_fixed


### Get predicted smooth effects (log scale) ---------------------------
pred_terms <- predict(m_rnd_ti, type = "terms")

# Look at the first few rows
head(pred_terms)
summary(pred_terms)

hist(pred_terms["s(prcp)",])

vis.gam(m_rnd_ti, view = c("prcp", "tmp"), plot.type = "contour")
# quick plotting ---------------------------
plot(m_rnd_ti_fixed, page = 1, shade = T)

m<- m_alt #m_rnd_ti         

clay_effect <- ggpredict(m, terms = "clay_extract")
dist_edge_effect <- ggpredict(m, terms = c("distance_edge"))
dist_sev_effect <- ggpredict(m, terms = c("disturbance_severity"))
dist_sev_effect_int <- ggpredict(m, terms = c('distance_edge', "disturbance_severity[0.3, 0.99]"))
tmp_prcp_effect <- ggpredict(m, terms = c("prcp", "tmp[8,10]"))
tmp_prcp_effect <- ggpredict(m, terms = c("drought_prcp", "tmp_z[0,1.5]"))

#manag_effect <- ggpredict(m, terms = c("protection_intensity", "country_pooled[CZ]"))


plot(clay_effect)
plot(dist_edge_effect)
plot(dist_sev_effect)
plot(dist_sev_effect_int)
plot(tmp_prcp_effect)
#plot(manag_effect)



## export final drivers model: only fixed effectd  -----------------------------------------------------------------

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
                  pred.labels = c("Intercept", pred_labels), # Replace smooth term labels
                  dv.labels = paste0("Explained Deviance: ", round(100 * summary(fin.m.reg.density)$dev.expl, 2), "%"), 
                  file = "outTable/full_drivers_reg_stem_density_fixed.doc")


#### with random effects ------------------

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
sjPlot::tab_model(m_rnd_ti,
                  show.re.var = TRUE,        # Show the variance components
                  pred.labels = c("Intercept", pred_labels), # Replace smooth term labels
                  dv.labels = paste0("Explained Deviance: ", round(100 * summary(m_rnd_ti)$dev.expl, 2), "%"), 
                  file = "outTable/full_drivers_reg_stem_density_full_random.doc")

## Disturbance severity: balanced design: --------------------------------
hist(df_stem_regeneration2$disturbance_severity)
hist(df_stem_regeneration2$distance_edge)

table(df_stem_regeneration2$disturb_sev_cl_simpl)[2]
table(df_stem_regeneration2$disturb_sev_cl3)
# min number of points in low category is 31
min_n = 175

# Randomly sample equal-sized subsets for each severity level
df_balanced2 <- df_stem_regeneration2 %>%
  group_by(disturb_sev_cl2) %>%
  sample_n(175) %>%   # from only 2 categories: low and high, low = 175
  ungroup()
# run for 3 caterories: low, medium, high
df_balanced3 <- df_stem_regeneration2 %>%
  group_by(disturb_sev_cl3) %>%
  sample_n(31) %>% # low category = 31
  ungroup()

# Fit the GAM model again on the balanced dataset
# Fit the GAM model again on the balanced dataset
m_balanced3 <- gam(
  stem_regeneration ~ s(prcp, k = 5) + s(tmp, k = 5) +
    te(disturbance_severity, distance_edge, k = 5) +
    ti(prcp, tmp, k = 5) +
    s(management_intensity, by = country_pooled, k = 4) +
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +
    s(x, y),
  family = tw(),
  method = "REML",
  data = df_balanced3
)



m_balanced2 <- gam(
  stem_regeneration ~ s(prcp, k = 5) + s(tmp, k = 5) +
    te(disturbance_severity, distance_edge, k = 5) +
    ti(prcp, tmp, k = 5) +
    s(management_intensity, by = country_pooled, k = 4) +
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +
    s(x, y),
  family = tw(),
  method = "REML",
  data = df_balanced2
)

# Compare the original and balanced models
summary(m_balanced3)
summary(m_balanced2)
summary(m_int_sev_edge_full_te_comb)  # Original model

hist(df_stem_regeneration2$management_intensity)

# Histogram for management intensity by country
ggplot(df_stem_regeneration2, aes(x = management_intensity)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = 
                   ) +
  facet_wrap(~ country_pooled, scales = "free_y") +
  labs(
    title = "Distribution of Management Intensity by Country",
    x = "Management Intensity",
    y = "Count"
  ) +
  theme_minimal()

ggplot(df_stem_regeneration2, aes(x = salvage_intensity)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = 
                   ) +
  facet_wrap(~ country_pooled, scales = "free_y") +
  labs(
    title = "Distribution of Salvage Intensity by Country",
    y = "Count"
  ) +
  theme_minimal()


ggplot(df_stem_regeneration2, aes(x = protection_intensity)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = 
                   ) +
  facet_wrap(~ country_pooled, scales = "free_y") +
  labs(
    title = "Distribution of Protection Intensity by Country",
    y = "Count"
  ) +
  theme_minimal()


# 6. interpret the results:  drivers -----------------------
# average increase in stem density per tmp, prcp and tehir interaction 
# Predict stem density at mean precipitation for different temperatures


## Predict effect of temperature (per 1°C increase)------------------
temp_pred <- ggpredict(fin.m.reg.density, terms = "tmp [7:12]")
temp_diff <- diff(temp_pred$predicted)

#Predict increase bfor TMP 
(avg_tmp_increase <- mean(temp_diff))
(avg_tmp <- mean(temp_pred$predicted))
(ci_tmp_low  <- mean(temp_pred$conf.low))
(ci_tmp_high <- mean(temp_pred$conf.high))
(avg_percent_increase <- (avg_tmp_increase / avg_tmp) * 100)

# Predict effect of precipitation (per 100mm increase)
prcp_pred <- ggpredict(fin.m.reg.density, terms = "prcp [500:1700 by=100]")
prcp_diff <- diff(prcp_pred$predicted)

(avg_prcp_increase <- mean(prcp_diff))
(avg_prcp          <- mean(prcp_pred$predicted))
(ci_prcp_low       <- mean(prcp_pred$conf.low))
(ci_prcp_high      <- mean(prcp_pred$conf.high))
(avg_percent_increase <- (avg_prcp_increase / avg_prcp) * 100)

# Predict interaction effect of temperature and precipitation
interaction_pred <- ggpredict(fin.m.reg.density, 
                              terms = c("prcp [500:1700 by=100]", "tmp [8,10]"))


# Separate predictions for each temperature
pred_8 <- interaction_pred %>% filter(group == "8")
pred_10 <- interaction_pred %>% filter(group == "10")

avg_prcp_level <- mean(pred_10$x) # same as avg_prcp_level <- mean(pred_8$x)


# Calculate the slope (change in density per 100mm increase in precipitation)
slope_8 <- diff(pred_8$predicted) / diff(pred_8$x)
slope_10 <- diff(pred_10$predicted) / diff(pred_10$x)

# Calculate the average slope for each temperature
avg_slope_8 <- mean(slope_8)
avg_slope_10 <- mean(slope_10)

# Mean stem density values for 8°C and 10°C
avg_stem_density_8 <- mean(pred_8$predicted)
avg_stem_density_10 <- mean(pred_10$predicted)

# Calculate the average percentage increase per 1°C
avg_interaction_percent_change <- ((avg_slope_10 - avg_slope_8) / avg_slope_8) * 100

# Display result
cat(sprintf(
  "On average, at 10°C, stem density increases %.2f%% faster with increasing precipitation compared to 8°C. ",
  avg_interaction_percent_change
))

cat(sprintf(
  "Mean stem density increased from %.0f stems/ha at 8°C to %.0f stems/ha at 10°C.\n",
  avg_stem_density_8, avg_stem_density_10
))


### distrubanec severity and distance edge -------------------------------------

# Predict effect of distance from edge
distance_pred <- ggpredict(fin.m.reg.density, terms = "distance_edge [50,250]")

# Extract predicted values
distance_50 <- distance_pred$predicted[distance_pred$x == 50]
distance_250 <- distance_pred$predicted[distance_pred$x == 250]

percent_drop <- ((distance_50 - distance_250) / distance_50) * 100
percent_drop

# Extract CIs
ci_50_low  <- distance_pred$conf.low[distance_pred$x == 50]
ci_50_high <- distance_pred$conf.high[distance_pred$x == 50]

ci_250_low  <- distance_pred$conf.low[distance_pred$x == 250]
ci_250_high <- distance_pred$conf.high[distance_pred$x == 250]

cat(sprintf(
  "At 50 m, predicted stem density was %.0f n ha⁻¹ (95%% CI: %.0f–%.0f), while at 250 m it was %.0f n ha⁻¹ (95%% CI: %.0f–%.0f).\n",
  distance_50, ci_50_low, ci_50_high,
  distance_250, ci_250_low, ci_250_high
))


## Predict stem density at disturbance severity levels 0.3 and 0.9 ---------------
severity_pred <- ggpredict(fin.m.reg.density, terms = "disturbance_severity [0.3,0.9]")

# Extract predicted values and confidence intervals
density_03     <- severity_pred$predicted[severity_pred$x == 0.3]
ci_03_low      <- severity_pred$conf.low[severity_pred$x == 0.3]
ci_03_high     <- severity_pred$conf.high[severity_pred$x == 0.3]

density_09     <- severity_pred$predicted[severity_pred$x == 0.9]
ci_09_low      <- severity_pred$conf.low[severity_pred$x == 0.9]
ci_09_high     <- severity_pred$conf.high[severity_pred$x == 0.9]

# Calculate percent increase
percent_increase <- ((density_09 - density_03) / density_03) * 100

# Print formatted sentence
cat(sprintf(
  "Stem density increased with disturbance severity.\n"
))
cat(sprintf(
  "At a severity of 0.3, predicted stem density was %.0f n ha⁻¹ (95%% CI: %.0f–%.0f),\n",
  density_03, ci_03_low, ci_03_high
))
cat(sprintf(
  "while at 0.9 it was %.0f n ha⁻¹ (95%% CI: %.0f–%.0f), representing an increase of %.0f%%.\n",
  density_09, ci_09_low, ci_09_high, percent_increase
))





### test for spatial autocorrelation: ------------------------------------------

# 

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


# 7.  Plot: Drivers  ---------------------------------------------------------------------------
y_lab = expression("Stem density [1000 n ha"^{-1}*"]")

summary(fin.m.reg.density)
anova.gam(fin.m.reg.density)

# show only 95% quatile for stem density
# Define the quantiles for stem_density and tmp columns
quantiles_stem_density99 <- quantile(df_stem_regeneration2$stem_regeneration  , 
                                    probs = c(0, 0.99), na.rm = TRUE)
# 
quantiles_stem_density90 <- quantile(df_stem_regeneration2$stem_regeneration  , 
                                     probs = c(0, 0.90), na.rm = TRUE)


# Filter the DataFrame to keep rows within these quantile ranges
filtered_df_plot99 <- df_stem_regeneration2 %>%
   dplyr::filter(stem_regeneration     >= quantiles_stem_density99[1] & stem_regeneration  <= quantiles_stem_density99[2])

filtered_df_plot90 <- df_stem_regeneration2 %>%
  dplyr::filter(stem_regeneration     >= quantiles_stem_density90[1] & stem_regeneration  <= 
                  quantiles_stem_density90[2])


# ,tmp >= quantiles_tmp[1] & tmp <= quantiles_tmp[2]

# Display the filtered data
#filtered_df_plot99
#summary(filtered_df_plot99$stem_regeneration)



## DRivers effect plorts make plot manually:  --------------------------

m <- fin.m.reg.density #m_int_sev_edge_full_te_comb # m_int_res_edge_full_te_comb #  m_int_res_edge_full_te
k.check(m)
summary(m)


# test 
# Generate predictions using ggpredict
# Interaction 1: Precipitation and Temperature
pred1 <- ggpredict(m, terms = c("prcp", "tmp [8,10]"))

pred1_df <- as.data.frame(pred1)
pred1_df$group <- as.numeric(as.character(pred1_df$group)) 
pred1_df$group <- factor(pred1_df$group)


# Interaction 2: Distance to Edge and Disturbance Severity
pred2 <- ggpredict(m, terms = c("distance_edge", "disturbance_severity [0.3,0.9]"))
# Example: Convert disturbance_severity to percent
pred2_df <- as.data.frame(pred2)
pred2_df$group <- as.numeric(as.character(pred2_df$group)) * 100
pred2_df$group <- factor(pred2_df$group)

my_colors_interaction <- c( "#FDAE61", "#A50026")
# Plot the first interaction
p1 <- 
  ggplot(pred1, aes(x = x, y = predicted/1000)) +
  geom_point(data = filtered_df_plot99, 
              aes( x = prcp, y = stem_regeneration/1000), 
              size = 1, alpha = 0.2, color = 'grey', pch = 16) +
  geom_line(linewidth = 1, aes(color = group) ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000, fill = group), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  scale_fill_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  theme_classic() +
  labs(x = "Precipitation [mm]", 
       y = "Regeneration stem density [#*1000/ha]", title = "p=0.005", 
      # linetype =  "Temperature [°C]"
       ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),       # Title size
    axis.title = element_text(size = 8),                   # Axis title size
    axis.text = element_text(size = 8),                    # Axis text size
    legend.key.size = unit(0.5, "cm"),                     # Legend key size
    legend.text = element_text(size = 8),                  # Legend text size
    legend.title = element_text(size = 8),                 # Legend title size
    legend.position = c(0.05, 0.9),
    legend.justification = c(0.05, 0.9)
  )

p1
# Plot the second interaction
p2 <- ggplot(pred2_df, aes(x = x, y = predicted/1000, color = group)) +
  geom_jitter(data = filtered_df_plot90, 
              aes( x = distance_edge, y = stem_regeneration/1000), 
              size = 1, alpha = 0.2, color = 'grey',
              width = 7,
              height = 1, pch = 16) +
  geom_line(linewidth = 1, aes(color = group) ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000, fill = group), alpha = 0.2, color = NA) +
  scale_color_manual(values = my_colors_interaction, 
                     name = "Disturbance\nseverity",
                     labels = c("Low", "High")) +
  scale_fill_manual(values = my_colors_interaction, 
                    name = "Disturbance\nseverity",
                    labels = c("Low", "High")) +
  theme_classic() +
  #ylim(0,20) +
  labs(x = "Distance to edge [m]", y = "", title = "p=0.751") +
  theme(
    axis.title = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 8),       # Title size
    axis.title.y = element_blank(),                   # Axis title size
    axis.text = element_text(size = 8),                    # Axis text size
    legend.key.size = unit(0.5, "cm"),                     # Legend key size
    legend.text = element_text(size = 8),                  # Legend text size
    legend.title = element_text(size = 8),                 # Legend title size
    legend.position = c(0.05, 0.9),
    legend.justification = c(0.05, 0.9)
  )
p2

p_combined_int <- ggarrange(p1, p2, 
                            labels = c("[a]","[b]"), 
                            align = 'hv',
                            font.label = list(size = 8, face = "plain")) # Specify plain font style)

(p_combined_int)
# Save the combined plot
ggsave('outFigs/fig_regen_int_drivers.png', plot = p_combined_int, 
       width = 6, height = 3.1, bg = 'white')


### p combined no points --------------------------

# Plot the first interaction
p1 <- 
  ggplot(pred1, aes(x = x, y = predicted/1000)) +
   geom_line(linewidth = 1, aes(color = group) ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000, fill = group), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  scale_fill_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  theme_classic() +
  labs(x = "Precipitation [mm]", 
       y = "Regeneration stem density [#*1000/ha]", title = "p=0.005", 
       # linetype =  "Temperature [°C]"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 8),       # Title size
    axis.title = element_text(size = 8),                   # Axis title size
    axis.text = element_text(size = 8),                    # Axis text size
    legend.key.size = unit(0.5, "cm"),                     # Legend key size
    legend.text = element_text(size = 8),                  # Legend text size
    legend.title = element_text(size = 8),                 # Legend title size
    legend.position = c(0.05, 0.9),
    legend.justification = c(0.05, 0.9)
  )

p1
# Plot the second interaction
p2 <- ggplot(pred2_df, aes(x = x, y = predicted/1000, color = group)) +
  geom_line(linewidth = 1, aes(color = group) ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000, fill = group), alpha = 0.2, color = NA) +
  scale_color_manual(values = my_colors_interaction, 
                     name = "Disturbance\nseverity",
                     labels = c("Low", "High")) +
  scale_fill_manual(values = my_colors_interaction, 
                    name = "Disturbance\nseverity",
                    labels = c("Low", "High")) +
  theme_classic() +
  #ylim(0,20) +
  labs(x = "Distance to edge [m]", y = "", title = "p=0.751") +
  theme(
    axis.title = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 8),       # Title size
    axis.title.y = element_blank(),                   # Axis title size
    axis.text = element_text(size = 8),                    # Axis text size
    legend.key.size = unit(0.5, "cm"),                     # Legend key size
    legend.text = element_text(size = 8),                  # Legend text size
    legend.title = element_text(size = 8),                 # Legend title size
    legend.position = c(0.05, 0.9),
    legend.justification = c(0.05, 0.9)
  )
p2

p_combined_int_no_points <- ggarrange(p1, p2, 
                            labels = c("[a]","[b]"), 
                            align = 'hv',
                            font.label = list(size = 8, face = "plain")) # Specify plain font style)

# Save the combined plot
ggsave('outFigs/fig_regen_int_drivers_no_points.png', plot = p_combined_int_no_points, 
       width = 6, height = 3.1, bg = 'white')

          
    

### Plot 4 drivers: clay, edge, severity, prcp_tmp ------------------------------------------
# Generate predictions using ggpredict
summary(m)

# Interaction 1: Precipitation and Temperature
pred_prcp       <- ggpredict(m, terms = c("prcp"))
pred_tmp        <- ggpredict(m, terms = c("tmp"))
pred_tmp_prcp   <- ggpredict(m, terms = c("prcp", "tmp [8,10]"))
pred_dist_edge  <- ggpredict(m, terms = c("distance_edge[50:250]"))
pred_dist_sever <- ggpredict(m, terms = c("disturbance_severity"))
pred_clay       <- ggpredict(m, terms = c("clay_extract"))

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

my_colors_interaction <- c("#FDAE61", "#A50026") 
my_color_main_effects <- "grey" # "#006837"    

# Plot precipitation
p.prcp <- ggplot(pred_prcp, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), alpha = 0.2,fill = my_color_main_effects) +
  labs(x = "Precipitation [mm]", 
       y = y_lab, 
       title = "p<0.001") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p.prcp

# Plot precipitation
p.tmp <- ggplot(pred_tmp, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), alpha = 0.2,fill = my_color_main_effects) +
  labs(x = "Temperature [°C]",  
       y = y_lab, 
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
       y =y_lab,  
       title = "p=0.006", 
       # linetype =  "Temperature [°C]"
  ) +
  my_theme_drivers

p1

# Plot distance to edge
p2 <- ggplot(pred_dist_edge, aes(x = x, y = predicted/1000)) +
   geom_line(linewidth = 1, color = my_color_main_effects) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), 
              alpha = 0.2,fill = my_color_main_effects) +
  labs(x = "Distance to edge [m]",  
       y = y_lab, 
       title = "p=0.023") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p2

# "#4D4D4D"
# "#666666"
# my_color_main_effects


# Plot the second interaction
p3 <- ggplot(pred_dist_sever, aes(x = x*100, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), fill =my_color_main_effects, alpha = 0.2, color = NA) +
  theme_classic() +
 # ylim(0,15) +
  labs(x = "Disturbance severity [%]", 
       y = "", title = "p=0.001") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p3


# PlotCaly
p4 <- ggplot(pred_clay, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000),  
              fill = my_color_main_effects, alpha = 0.2, color = NA) +
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
ggsave('outFigs/fig_regen_int_drivers_no_points.png', plot = p_combined_int_no_points, 
       width = 5.5, height = 5.5, bg = 'white')



## all predictors for supplemeny  ------------------------

p_combined_int_no_points_supplem <- ggarrange(p.prcp, p.tmp, p1, p4, p2, p3,
                                      labels = c("[a]","[b]", "[c]","[d]", "[e]","[f]"), 
                                      ncol = 3, nrow = 2,
                                      align = 'hv',
                                      font.label = list(size = 8, face = "plain")) # Specify plain font style)

p_combined_int_no_points_supplem

# Save the combined plot
ggsave('outFigs/fig_p_combined_int_no_points_supplem.png', plot = p_combined_int_no_points_supplem, 
       width = 8, height = 5.5, bg = 'white')



## 6. Wilcox: Boxplot sites differences: delayed vs advaced ----------------------------

# two categories: count how many plots i have?

#prop.table(table(df_fin$adv_delayed_wider))
prop.table(table(df_fin$adv_delayed))
table(df_fin$adv_delayed)




# get indicators for delayed/advanced: ------------------------------------------
df_delayed_advanced <- df_fin %>% 
  dplyr::select(site, country, adv_delayed)


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
                sd_grw_anm_prcp,
                sd_grw_anm_tmp,
                tmp_z,
                prcp_z,
                #max_grw_anm_prcp,
               # max_grw_anm_tmp,
                av.nitro, 
               depth_extract,
               slope, 
               elevation,
               sand_extract,
               clay_extract,
               
               #management_intensity,
               #salvage_intensity,
               #protection_intensity,
              # sapl_juv_ratio,
                #cv_t2m, 
               # cv_tp,
                
                #spei12,
                #spei1,
                #spei6,
                drought_spei1,
                drought_spei3,
                drought_spei12,
                #drought_tmp,
                drought_prcp,
                disturbance_severity,
                distance_edge, 
              #  time_since_disturbance,
                adv_delayed) %>% 
  #dplyr::select(all_of(variables_to_plot), adv_delayed) %>% 
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

print(summary_stats_narrow_sub)


# list groups to pairwise comparison
comparisons_3grps <- list(c("Delayed", "Other"), 
                          c("Delayed", "Advanced"), 
                          c("Other", "Advanced"))

comparisons_2grps <- list(c("Delayed", "Advanced"))


# # Compute Wilcoxon effect sizes for each variable
# library(rstatix)
# 
# effect_sizes <- df_long_narrow %>%
#   group_by(Variable) %>%
#   wilcox_effsize(Value ~ adv_delayed, comparisons = comparisons)
# 
# # most of my effect sizes (differences between two groups) are rather small, the highest effect has drought_prcp
# # boot strap analysis: 
# 
# library(boot)
# boot_effsize <- boot(data = df_long_narrow, statistic = function(data, indices) {
#   sample_data <- data[indices, ]
#   wilcox_effsize(sample_data$Value ~ sample_data$adv_delayed)
# }, R = 1000)

# Calculate mean ± SD and median ± IQR for each Variable per group - all variables
summary_stats_adv_delayed <- df_long_narrow %>%
  group_by(Variable, adv_delayed) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    IQR = IQR(Value, na.rm = TRUE)#,
    #Min = min(Value, na.rm = TRUE),
    #Max = max(Value, na.rm = TRUE),
    #n = n()  # Number of observations per group
  ) %>%
  ungroup()


print(summary_stats_adv_delayed, n = 40)


# Plot using ggboxplot
p_boxplot_wilcox_narrow <- 
  ggboxplot(df_long_narrow, x = "adv_delayed", y = "Value", 
                              fill = "adv_delayed", 
                              palette = c("lightblue", "red", "green"),
                              facet.by = "Variable", scales = "free_y", 
                              ylab = "Values", xlab = "Regeneration Status",
                              #outlier.size = .2,
                              outlier.shape = NA,
                              size = 0.2) +
  labs(title = "",
       x = "Reg. Status", y = "Vals")+
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
    panel.border = element_rect(colour = 
                                  , fill = NA, linewidth = 0.5)  # Add a square border around the plots
  )  +
  stat_compare_means(comparisons = comparisons_3grps, method = "wilcox.test", 
                     label = "p.signif", #"p.format",  
                     size = 2, label.x = 1.5) +
  # Add mean dots
    geom_point(data = summary_stats_adv_delayed, 
               aes(x = adv_delayed, y = Mean, group = Variable), 
               shape = 21, 
               fill = "white",        # Added missing fill value
               color = "black",       # Added missing color value
               size = 1.5, 
               inherit.aes = FALSE)

p_boxplot_wilcox_narrow

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


## Supplement Wilcox: ---make a final plot: only prcp,  tmp_z, drought_spei1, - 3  groups -------------------
# manually 

# Example: Create individual plots
p.prcp <-   df_long_narrow_sub %>% 
  dplyr::filter(Variable == "prcp") %>% 
  ggboxplot(x = "adv_delayed", 
            y = "Value", 
            fill = "adv_delayed", 
            palette = c("#A50026", "#FDAE61", "#006837"),
            ylab = "Precipitation [mm]", 
            xlab = "",
            outlier.shape = NA,
            #outlier.size = .2,
            size = 0.2) +
  stat_compare_means(comparisons = comparisons_3grps, 
                     method = "wilcox.test", 
                     label =  'p.format', #"p.signif",  
                     size = 3,
                     label.y = c(1500, 
                                 1400, 
                                 1300)) +  # Standard Wilcoxon test comparison lines
  my_theme() +
  # Add mean dots
  geom_point(data =  subset(summary_stats_narrow_sub, Variable == "prcp"), 
             aes(x = adv_delayed, y = Mean, group = adv_delayed), 
             shape = 21, fill = "red", color = "red", size = 1.5, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(300,1600))  #

p.prcp


p.tmp <- 

    df_long_narrow_sub %>% 
  dplyr::filter(Variable == "tmp") %>% 
    ggboxplot(x = "adv_delayed", 
              y = "Value", 
                   fill = "adv_delayed", 
                   palette = c("#A50026", "#FDAE61", "#006837"),
                   ylab = "Temperature [°C]",
                   xlab = "",
                   outlier.shape = NA,
                   #outlier.size = .2,
                   size = 0.2) +
      stat_compare_means(comparisons = comparisons_3grps,
                         method = "wilcox.test",
                         label =  'p.format', #"p.signif",
                         size = 3,
                         label.y = c(13.4,
                                     12.8,
                                     12.2)) +
  my_theme() +
      # Add mean dots
      geom_point(data =  subset(summary_stats_narrow_sub, Variable == "tmp"), 
                 aes(x = adv_delayed, y = Mean, group = adv_delayed), 
                 shape = 21, fill = "red", color = "red", size = 1.5, inherit.aes = FALSE) +
    coord_cartesian(ylim = c(7,14)) 
  
p.tmp

p.spei1 <-
  df_long_narrow_sub %>% 
  dplyr::filter(Variable == "drought_spei1") %>% 
  ggboxplot(x = "adv_delayed", 
            y = "Value", 
            fill = "adv_delayed", 
            palette = c("#A50026", "#FDAE61", "#006837"),
            ylab = "SPEI-1 Drought Index [dim.]", 
            xlab = "",
            outlier.shape = NA,
            #outlier.size = .2,
            size = 0.2) +
  stat_compare_means(comparisons = comparisons_3grps, 
                     method = "wilcox.test", 
                     label =  'p.format', #"p.signif",  
                     size = 3,
                     label.y = c(-0.7, 
                                 -0.75, 
                                 -0.8)) +  # Standard Wilcoxon test comparison lines
  my_theme() +
  # Add mean dots
  geom_point(data =  subset(summary_stats_narrow_sub, Variable == "drought_spei1"), 
             aes(x = adv_delayed, y = Mean, group = adv_delayed), 
             shape = 21, fill = "red", color = "red", size = 1.5, inherit.aes = FALSE) +

  coord_cartesian(ylim = c(-1.15,-0.65)) 

p.spei1




# Example: Create individual plots
p.clay <-   
  df_long_narrow_sub %>% 
  dplyr::filter(Variable == "clay_extract") %>% 
  ggboxplot(x = "adv_delayed", 
            y = "Value", 
            fill = "adv_delayed", 
            palette = c("#A50026", "#FDAE61", "#006837"),
            ylab = "Clay [%]", 
            xlab = "",
            outlier.shape = NA,
            #outlier.size = .2,
            size = 0.2) +
  stat_compare_means(comparisons = comparisons_3grps, 
                     method = "wilcox.test", 
                     label =  'p.format', #"p.signif",  
                     size = 3,
                     label.y = c(42, 
                                 38, 
                                 35)) +  # Standard Wilcoxon test comparison lines
  my_theme() +
  # Add mean dots
  geom_point(data =  subset(summary_stats_narrow_sub, Variable == "clay_extract"), 
             aes(x = adv_delayed, y = Mean, group = adv_delayed), 
             shape = 21, fill = "red", color = "red", size = 1.5, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0,45))  #

p.clay

# Example: Create individual plots
p.disturbance <-   df_long_narrow_sub %>% 
  dplyr::filter(Variable == "disturbance_severity") %>% 
  ggboxplot(x = "adv_delayed", 
            y = "Value", 
            fill = "adv_delayed", 
            palette = c("#A50026", "#FDAE61", "#006837"),
            ylab = "Disturbance severity [dim.]", 
            xlab = "",
            outlier.shape = NA,
            #outlier.size = .2,
            size = 0.2) +
  stat_compare_means(comparisons = comparisons_3grps, 
                     method = "wilcox.test", 
                     label =  'p.format', #"p.signif",  
                     size = 3,
                     label.y = c(1.2, 
                                 1.1, 
                                 1)) +  # Standard Wilcoxon test comparison lines
  my_theme() +
  # Add mean dots
  geom_point(data =  subset(summary_stats_narrow_sub, Variable == "disturbance_severity"), 
             aes(x = adv_delayed, y = Mean, group = adv_delayed), 
             shape = 21, fill = "red", color = "red", size = 1.5, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0,1.3))  #

p.disturbance


# Example: Create individual plots
p.distance_edge <-   
  df_long_narrow_sub %>% 
  dplyr::filter(Variable == "distance_edge") %>% 
  ggboxplot(x = "adv_delayed", 
            y = "Value", 
            fill = "adv_delayed", 
            palette = c("#A50026", "#FDAE61", "#006837"),
            ylab = "Distance to edge [m]", 
            xlab = "",
            outlier.shape = NA,
            #outlier.size = .2,
            size = 0.2) +
  stat_compare_means(comparisons = comparisons_3grps, 
                     method = "wilcox.test", 
                     label =  'p.format', #"p.signif",  
                     size = 3,
                     label.y = c(140, 
                                 130, 
                                 120)) +  # Standard Wilcoxon test comparison lines
  my_theme() +
  # Add mean dots
  geom_point(data =  subset(summary_stats_narrow_sub, Variable == "distance_edge"), 
             aes(x = adv_delayed, y = Mean, group = adv_delayed), 
             shape = 21, fill = "red", color = "red", size = 1.5, inherit.aes = FALSE) +
  coord_cartesian(ylim = c(0,150))  #

p.distance_edge



# Combine all plots into a single figure
wilcox_plot_out_clay <- ggarrange(p.prcp, p.tmp, p.spei1,p.clay,
                        ncol = 2, 
                        nrow = 2, 
                        labels = c("[a]", "[b]","[c]","[d]"), font.label = list(size = 8, face = "plain"))

# Display the final merged plot
wilcox_plot_out_clay


# Save the plot as an SVG file
ggsave(filename = "outFigs/wilcox_plot_out_clay.png", plot = wilcox_plot_out_clay, 
       device = "png", width = 5, height = 5.5, dpi = 300, bg = 'white')

# Save the plot as an SVG file
ggsave(filename = "outFigs/wilcox_plot_out_clay.svg", plot = wilcox_plot_out_clay, 
       device = "svg", width = 5, height = 5.5, dpi = 300, bg = 'white')



# Wilcox supplement ------------------
# alternative Wilcopx with protection intensity
# Combine all plots into a single figure
wilcox_plot_supplement<- ggarrange(p.prcp, p.tmp, p.spei1, p.clay, #p.manag, p.protection,
                                   p.disturbance, p.distance_edge,
                              ncol = 3, nrow = 2, labels = c("[a]", "[b]","[c]","[d]", '[e]','[f]'), #'[g]','[h]'
                              font.label = list(size = 8, face = "plain"))

# Save the plot as an SVG file
ggsave(filename = "outFigs/wilcox_vars_supplement.png", plot = wilcox_plot_supplement, 
       device = "png", width = 5.5, height = 4, dpi = 300, bg = 'white')



### FINAL wilcox: remove 'other' category -----
# filter the tables
# !!! START
summary_stats_narrow_sub2 <- summary_stats_narrow_sub %>%
  dplyr::filter(adv_delayed != "Other") %>%
  droplevels()

df_long_narrow_sub2 <- df_long_narrow_sub %>%
  dplyr::filter(adv_delayed != "Other") %>%
  droplevels()


# > species_colors
# piab      fasy      quro      pisy      soau      acps      potr      abal 
# "#006837" "#229C52" "#74C364" "#B7E075" "#E9F6A2" "#FEEDA2" "#FDBE6E" "#F67B49" 
# besp      lade 
 #"#DA362A" "#A50026" 
# 
# make a function  # "#FDAE61", 
plot_variable_box <- function(var_name, 
                              y_label, 
                              y_limits, 
                              label_y_values, 
                              palette = c(  "#F67B49", "#74C364")) {
  
  data_filtered <- df_long_narrow_sub2 %>%
    dplyr::filter(Variable == var_name)
  
  # Base plot
  p <- ggboxplot(data_filtered,
                 x = "adv_delayed", 
                 y = "Value", 
                 fill = "adv_delayed", 
                 palette = palette,
                 ylab = y_label, 
                 xlab = "",
                 outlier.shape = NA,
                 size = 0.2) +
    
    # Jittered raw data points
    # geom_jitter(data = data_filtered,
    #             aes(x = adv_delayed, y = Value),
    #            color = 'gray30',
    #              width = 0.2,
    #             alpha = 0.4,
    #             size = 0.8,
    #             inherit.aes = FALSE) +
    # 
    # Wilcoxon comparison
    stat_compare_means(comparisons = comparisons_2grps, 
                       method = "wilcox.test", 
                       label = 'p.format',
                       size = 3,
                       label.y = label_y_values) +
    
    # Theme
    my_theme() +
    
    # Red mean points
    geom_point(data = subset(summary_stats_narrow_sub2, Variable == var_name), 
               aes(x = adv_delayed, y = Mean, group = adv_delayed), 
               shape = 21, fill = "black", color = "black", size = 1.5, inherit.aes = FALSE) +
    
    # Zoomed y-axis
    coord_cartesian(ylim = y_limits)
  
  return(p)
}

p.spei1.fin <- plot_variable_box(
  var_name = "drought_spei1",
  y_label = "SPEI-1 [dim.]",
  y_limits = c(-1.2, -0.6),
  label_y_values = c( -0.62, -0.6)
)

p.prcp.fin <- plot_variable_box(
  var_name = "prcp",
  y_label = "Precipitation [mm]",
  y_limits = c(300, 1600),
  label_y_values = c(1500, 1400, 1300)
)

p.clay.fin <- plot_variable_box(
  var_name = "clay_extract",
  y_label = "Clay [%]",
  y_limits = c(0, 45),
  label_y_values = c(42, 38, 35)
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
ggsave(filename = "outFigs/wilcox_fin.png", plot = wilcox_fin, 
       device = "png", width = 6, height = 3, dpi = 300, bg = 'white')

# Save the plot as an SVG file
ggsave(filename = "outFigs/wilcox_fin.svg", plot = wilcox_fin, 
       device = "svg", width = 6, height = 3, dpi = 300, bg = 'white')


# use bayes statistics: more robust for uneven sample sizes ----------------------------
library("BayesFactor")
# Compute Bayesian ANOVA for each variable
bayesian_results <- df_long_narrow_sub2 %>%
  group_by(Variable) %>%
  summarise(
    BF = extractBF(anovaBF(Value ~ adv_delayed, data = cur_data()))$bf[1]  # Extract Bayes Factor
  ) %>% 
  arrange(decs(BF))


sjPlot::tab_df(bayesian_results,
               show.rownames = FALSE,
               file="outTable/bayesian_results_groups_diff.doc",
               digits = 1) 


# test if i got anova right --
df_sub <- df_long_narrow %>% 
  dplyr::filter(Variable == 'drought_prcp')

anovaBF(Value ~ adv_delayed, data = df_sub)



# Test 

df_fin %>% 
  group_by(adv_delayed) %>% 
  summarize(
    min_distance = min(distance_edge, na.rm = TRUE),
    max_distance = max(distance_edge, na.rm = TRUE),
    mean_distance = mean(distance_edge, na.rm = TRUE),
    median_distance = median(distance_edge, na.rm = TRUE),
    IQR_distance = IQR(distance_edge, na.rm = TRUE),
    sd_distance = sd(distance_edge, na.rm = TRUE),
    n = n()
  )

df_fin %>% 
  #group_by(adv_delayed) %>% 
  summarize(
    min_distance = min(distance_edge, na.rm = TRUE),
    max_distance = max(distance_edge, na.rm = TRUE),
    mean_distance = mean(distance_edge, na.rm = TRUE),
    median_distance = median(distance_edge, na.rm = TRUE),
    IQR_distance = IQR(distance_edge, na.rm = TRUE),
    sd_distance = sd(distance_edge, na.rm = TRUE),
    n = n()
  )


df_fin %>% 
  #group_by(adv_delayed) %>% 
  summarize(
    min_severity = min(disturbance_severity, na.rm = TRUE),
    max_severity = max(disturbance_severity, na.rm = TRUE),
    mean_severity = mean(disturbance_severity, na.rm = TRUE),
    median_severity = median(disturbance_severity, na.rm = TRUE),
    IQR_severity = IQR(disturbance_severity, na.rm = TRUE),
    sd_severity = sd(disturbance_severity, na.rm = TRUE),
    n = n()
  )





# 7. Analyze at country level: Country effect ----------------------------------------------------------
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

p.country.richness

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


### Get partial R2: define teh most important variable --------------------
# Perform an ANOVA to assess each term's contribution to the model

# Run ANOVA on GAM model
anova_results <- anova(fin.m.reg.density)

# Convert the ANOVA results into a dataframe
anova_df <- as.data.frame(anova_results)

# Add term names for better readability
anova_df <- tibble::rownames_to_column(anova_df, "Predictor")

# Rename columns for better interpretation
colnames(anova_df) <- c("Predictor", "EDF", "Ref. DF", "F-Value", "p-Value")

# Export as an APA-styled table using sjPlot
sjPlot::tab_df(anova_results,
               title = "outTable/ANOVA Results for GAM Model",
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = FALSE,
               file="outTable/anova_GAM.doc",
               digits = 3) 


## Save models --------------------------------------------------

# Save the model object and input data
save(fin.m.reg.density, df_stem_regeneration2,
     file = "outData/stem_density_models.RData")



## Save data ------------------------------------------------
fwrite(df_delayed_advanced, 'outTable/df_delayed_advanced.csv')
fwrite(df_stem_regeneration2, 'outTable/df_stem_regeneration2.csv')


# 5. Species climate suitability  -------------------------------------------

# overall richness from field data?

stem_dens_species_long_cluster %>% 
  dplyr::filter(stem_density > 0) %>% 
  distinct(Species) %>% 
  count()  #@ 35

# compare current species composition with future ones form wessely
# Need to adjust tree species naming a bit: 
#   we have less species then wessely: 
  # - eg we have besp, which can be betula pendula, betula pubescense
  # - same for populus (3), quercus (6), salix ...

head(stem_dens_species_long_cluster)

# get species merging table - acronyms from iLand
species_look_up_simple<- read.csv("rawData/tree_sp_simple.csv", sep = ';')

# get species merging table:   - merged my naming with Wessely species naming
species_look_up_full<- read.csv("rawData/tree_sp_field_wessely_merged.csv", sep = ';')


# identify what species are present per plot - use simple look up table, 
# add only latin names
present_species <- 
  stem_dens_species_long_cluster %>% 
  dplyr::filter(VegType != "Mature") %>% 
  ungroup() %>% 
  group_by(cluster, Species) %>% 
  summarize(sum_stem_density = sum(stem_density, na.rm = T )) %>% 
  mutate(presence = if_else(sum_stem_density > 0, 1, 0)) %>% 
 # dplyr::select(-sum_stem_density ) %>% 
  left_join(species_look_up_simple, by = c("Species" = "acc")) %>% 
  dplyr::select(-latin) 

present_species <- present_species %>% 
  dplyr::rename(acc = Species) %>% 
  dplyr::rename(site = cluster)

# read table with Wessely species - from species present/absent of clusters locations
# only species present on every perios (occurence = 8) are climatically suitable across whole century
future_species      <- fread('outTable/species_presence_clim_change.csv')
# future_species_full <- fread('outTable/species_presence_clim_change_full.csv')

unique(future_species$species )


# add acronyms and consider presence across several species: eg betula sp.
future_species_sum <- 
  future_species %>%
  # add naming
  left_join(species_look_up_full, by = c('species' = 'wessely')) %>%
 #   View()
  # group by species to allow occurence of species that have specified genus: eg betula sp.
  group_by(site, scenario, acc) %>% 
  # Summarize and set sum_presence to 1 if the sum is greater than 1 - 
  # this account for the fact that wessely can have two betulas, I have only 1
  summarize(sum_presence = pmin(sum(overall_presence), 1), .groups = 'drop')
  
wide_future_species <- future_species_sum %>%
  pivot_wider(names_from = scenario, values_from = sum_presence) 

# merge both tables: the presently recorded species and species under climate scenarios
df_compare_future_species <- wide_future_species %>% 
  left_join(present_species) %>%   # use left join to explude species that are recorded in field, but not present in Wessely database
  dplyr::rename(current = presence)

anyNA(df_compare_future_species)
length(unique(df_compare_future_species$acc))  # final length is 30 species: corss between wessely and observed field database

# how many species are present in my species database, but not in Wessely? 
species_out_of_wessely <- df_compare_future_species %>%
  dplyr::filter(is.na(rcp26)) %>%
  distinct(acc)  # Optional, if you want all unique rows with NA acc

# add country indication
df_compare_future_species <- df_compare_future_species %>%
  # Extract the first two characters of 'site' as 'region' and convert to integer
  mutate(region = as.integer(substr(site, 1, 2))) %>%
  # Left join with unique_regions_per_country to get country indication
  left_join(unique_regions_per_country, by = c("region" = "unique_regions"))

# evaluate presence vs absence of species per cluster
df_compare_future_species <- df_compare_future_species %>%
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


# Share of stems that is not suitable under CC scenarios?
# Convert suitability columns to logical: TRUE if unsuitable ("0"), FALSE otherwise

anyNA(df_compare_future_species)  

# get total sum of stems: acdjust for fact that ionly 30 species overlaps between 2 dababases
total_stems         <- sum(stem_dens_species_long_cluster$stem_density)
total_stems_wessely <- sum(df_compare_future_species$sum_stem_density )
total_stems_wessely
#6172500 - all stems, 6008250 over crossed databases 

# now, it is only climate suitability for regeneration!
df_compare_future_species_long_full <- df_compare_future_species %>%
  dplyr::select(site, 
                acc,
                sum_stem_density , 
                country_pooled, 
                suitability_rcp26, 
                suitability_rcp45, 
                suitability_rcp85) %>%
  pivot_longer(cols = starts_with("suitability_rcp"), 
               names_to = "scenario", 
               values_to = "suitability") %>%
  mutate(scenario = gsub("suitability_", "", scenario)) #


df_compare_future_species_long_full %>%  # Remove "suitability_" prefix
  group_by(suitability, scenario) %>% 
  dplyr::filter(suitability == "not_suitable") %>%  
  summarize(sum = sum(sum_stem_density)) %>% 
    mutate(share = round(sum/total_stems*100,1))
 
  # suitability  scenario     sum share
  # <chr>        <chr>      <dbl> <dbl>
  # 1 not_suitable rcp26    3838375  62.2
  # 2 not_suitable rcp45    4505375  73  
  # 3 not_suitable rcp85    5172625  83.8


# analysis over species: --------------------------
# how many plots of piab swill not be climatically suitable?


### across all species ---------------
# Filter for species with stem density > 0
species_suitability_summary <- df_compare_future_species_long_full %>%
  dplyr::filter(sum_stem_density > 0) %>%
  
  # Count total number of unique sites per species
  group_by(acc) %>%
  mutate(total_sites = n_distinct(site)) %>%
  
  # Filter to only suitable cases
  dplyr::filter(suitability == "suitable") %>%
  
  # Count number of suitable sites per species and scenario
  distinct(site, scenario, acc, total_sites) %>%
  count(acc, scenario, name = "n_suitable_sites") %>%
  
  # Join back total site counts
  left_join(
    df_compare_future_species_long_full %>%
      dplyr::filter(sum_stem_density > 0) %>%
      group_by(acc) %>%
      summarise(n_sites_with_species = n_distinct(site), .groups = "drop"),
    by = "acc"
  ) %>%
  
  # Calculate share
  mutate(n_unsuitable_sites = n_sites_with_species - n_suitable_sites,
         share_suitable     = round(n_suitable_sites / n_sites_with_species * 100,1),
         share_unsuitable   = round(n_unsuitable_sites / n_sites_with_species * 100,1),
         share_overall      = round(n_sites_with_species / total_sites * 100,1)) %>%
  dplyr::select(scenario, 
                share_overall,
                n_suitable_sites, 
                n_unsuitable_sites,     
                n_sites_with_species,   # number of sites that species occured on
                share_suitable,         
                share_unsuitable,
                acc)

# View result
species_suitability_summary

# identify the most affected species (highest share of unsuitable on plot in RCP45)
species_suitability_summary %>%
  dplyr::filter(scenario == "rcp45") %>%
  #dplyr::filter(acc %in% top_species_site_share) %>% 
  mutate(share_unsuitable = n_unsuitable_sites / n_sites_with_species * 100) %>%
  arrange(desc(share_overall)) %>% # , share_unsuitable
View()

# the most affected species by cliamte change
share_unsuitable_range <- species_suitability_summary %>%
  group_by(acc) %>%
  summarise(
    min_unsuitable = min(share_unsuitable),
    max_unsuitable = max(share_unsuitable),
    range_unsuitable = max_unsuitable - min_unsuitable,
    .groups = "drop"
  ) %>%
  arrange(desc(range_unsuitable))

share_unsuitable_range


## Plot level : Get only counts per country and plot, no need for specific species : ------------
# Convert to long format
df_suitability_long <- 
  df_compare_future_species %>%
  dplyr::select(site, 
                country_pooled, 
                suitability_rcp26, 
                suitability_rcp45, 
                suitability_rcp85) %>%
  pivot_longer(cols = starts_with("suitability_rcp"), 
               names_to = "scenario", 
               values_to = "suitability") %>%
  mutate(scenario = gsub("suitability_", "", scenario)) %>%  # Remove "suitability_" prefix
  group_by(site, country_pooled, suitability, scenario) %>% 
  summarize(freq = n()) %>% 
  dplyr::filter(suitability != 0) %>% 
  ungroup()

# how many plots overall do not have any tree species currently present  climatically suitable? 
df_suitability_long %>%
    group_by(site, scenario) %>%
    summarize(has_suitable = any(suitability == "suitable"), .groups = "drop") %>%
    filter(!has_suitable) %>%
    count(scenario, name = "n_clusters_without_suitable") %>% 
    mutate(share = n_clusters_without_suitable/849)
  
# scenario n_clusters_without_suitable share
# <chr>                          <int> <dbl>
#   1 rcp26                            327 0.385
# 2 rcp45                            433 0.510
# 3 rcp85                            524 0.617 
  
df_suitability_summary_plot <- df_suitability_long %>% 
  group_by( country_pooled, suitability, scenario) %>% 
  dplyr::summarize(mean = mean(freq, na.rm = T),
                   sd = sd(freq, na.rm = T)) %>% 
  mutate(cv = (sd / mean) * 100) %>% 
  mutate(suitability = factor(suitability, levels = c("not_suitable", "suitable", "novel"))) %>% 
  mutate(mean = ifelse(suitability == "not_suitable", -mean, mean)) %>% # change values to neagative for not_suitable species
  ungroup(.)


#df_suitability_summary_plot <- df_suitability_summary_plot[order(levels(df_suitability_summary_plot$suitability)),]
# make a barplot like this :

# figure out teh error bars!!! 
df_suitability_summary_plot2 <- df_suitability_summary_plot %>%
  arrange(country_pooled, scenario) %>%
  group_by(country_pooled, scenario) %>%
  mutate(
    # Compute cumulative positions only for suitable and novel
    cumulative_mean = ifelse(suitability == "suitable", cumsum(mean), mean),
    cumulative_mean = ifelse(suitability == "not_suitable", mean, -mean ),
    cumulative_mean = ifelse(suitability == "novel", mean, mean ),
    
    # Adjust error bars: Flip for "not_suitable"
    cumulative_upper = ifelse(suitability == "not_suitable", mean - sd, cumulative_mean + sd),
    cumulative_lower = ifelse(suitability == "not_suitable", mean + sd, cumulative_mean - sd)
  ) %>% 
  mutate(
    pos = cumsum(mean),
    upper = pos + sd/2,
    lower = pos - sd/2
  ) %>%
  ungroup()# with error bars: 

# add specific factor for x:
df_suitability_summary_plot2 <- df_suitability_summary_plot2 %>% 
  dplyr::filter(suitability !='novel') %>% 
  mutate(x_comb = factor(paste(country_pooled, scenario )))


# Test start -----!!!

# Create scenario labels for x-axis (only show for unique scenario values)
scenario_labels <- df_suitability_summary_plot2 %>%
  group_by(x_comb) %>%
  summarize(label = first(scenario)) %>%
  pull(label)

p.species.suitability_country <- 
  ggplot(df_suitability_summary_plot2, aes(x = scenario, 
                                           y = mean, 
                                           fill =  factor(suitability, 
                                           levels = c("not_suitable", "suitable", "novel")))) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  #  geom_errorbar(aes(ymin = lower, ymax = upper), 
  #          width = 0.2, color = "grey10") +  # Add error bars at midpoints
  theme_classic2() +
    scale_x_discrete(labels = scenario_labels) +
    facet_grid(. ~ country_pooled) +  # Facet by country, labels below
  labs(title = "",
       x = "",
       y = "Number of species [#]",
       fill = "Climate suitability") +
    scale_fill_manual(
      values = c("not_suitable" = "#FDBE6E", 
                 "suitable" = "#006837", 
                 "novel" = "#FEEDA2"),  # Custom colors
      labels = c("not_suitable" = "Not suitable", 
                 "suitable" = "Suitable", 
                 "novel" = "Novel")  # Custom legend labels
    )+
 # geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  theme(
    legend.position = 'right',
    legend.text = element_text(size = 8),  
    legend.title = element_text(size = 8),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),  # Rotate x-axis labels
    axis.text.y = element_text(size = 8),  # Keep y-axis labels non-italic
    axis.title = element_text(size = 8),    # Adjust axis title size
    strip.background = element_blank(),  # Removes the box around facet labels
    strip.text = element_text(face = "bold"),  # Keeps facet text bold
    strip.placement = "outside",  # Places facet labels outside the plot area
    panel.spacing = unit(0.3, "lines")  # Increases space between facet panels (countries)
  )  
p.species.suitability_country

#"#006837" "#229C52" "#74C364" "#B7E075" "#E9F6A2" "#FEEDA2" "#FDBE6E" "#F67B49" "#DA362A" "#A50026" 
## Country level: identify species suitability on country level: calculate richness per country at current, at RCP45 and compare both ------

df_suitab_sub45 <-df_compare_future_species %>% 
  dplyr::select(acc, suitability_rcp45, country_pooled)


df_species_current <- df_compare_future_species %>%
  #dplyr::filter(country_pooled == "AT") %>% 
  dplyr::filter(current == 1) %>% 
  dplyr::select(acc, country_pooled) %>% 
  distinct() %>% 
  mutate(current = 1)

df_species_rcp45 <- df_compare_future_species %>%
  #dplyr::filter(country_pooled == "AT") %>% 
  dplyr::filter(rcp45 == 1) %>% 
  dplyr::select(acc, country_pooled) %>% 
  distinct() %>% 
  mutate(rcp45 = 1)

df_species_comparison <- full_join(df_species_current, df_species_rcp45, 
                                   by = c("acc", "country_pooled")) %>%
  mutate(
    current = ifelse(is.na(current), 0, current),  # Fill NA with 0
    rcp45 = ifelse(is.na(rcp45), 0, rcp45)        # Fill NA with 0
  ) %>%
  mutate(  suitability_rcp45 = case_when(
    current == 1 & rcp45 == 1 ~ "suitable",
    current == 1 & rcp45 == 0 ~ "not_suitable",
    current == 0 & rcp45 == 1 ~ "novel",
    current == 0 & rcp45 == 0 ~ "0"
  ))

# summarize information fo species on country level: species counts
# Count occurrences of "suitable" and "not_suitable" per country
df_suitability_summary <- df_species_comparison %>%
  group_by(country_pooled, suitability_rcp45) %>%
  summarise(species_count = n(), .groups = "drop") #%>%
  #tidyr::pivot_wider(names_from = suitability_rcp45, values_from = species_count, values_fill = 0)


df_richness_country_simpl <- df_richness_country_current %>% 
  dplyr::select(-acc) %>% 
  distinct() %>% 
  rename(current_richness = richness)

# Rename columns for clarity
df_suitability_summary <- df_suitability_summary %>%
  left_join(df_richness_country_simpl) #%>% 
  #dplyr::filter(suitability_rcp45 != 'novel')


# Display result
print(df_suitability_summary)

ggplot(df_suitability_summary, aes(x = country_pooled, y = species_count, fill = suitability_rcp45)) +
  geom_bar(stat = "identity", position = 'fill') +  # Stacked bar plot
  theme_classic2() +
  labs(title = "Species Suitability per Country under RCP45",
       x = "",
       y = "Share of species [%]",
       fill = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  scale_fill_manual(values = c("suitable" = "#006837", "not_suitable" = "#FDBE6E", "novel" = "#FEEDA2"))  # Custom colors


# "#006837" "#229C52" "#74C364" "#B7E075" "#E9F6A2" "#FEEDA2" "#FDBE6E" "#F67B49" "#DA362A" "#A50026" 
# get just list of unique species per country under scenarios ---------------







# Convert to long format
df_long_suitability <- df_compare_future_species %>%
  pivot_longer(cols = starts_with("suitability_rcp"), 
               names_to = "scenario", 
               values_to = "suitability") %>%
  mutate(scenario = dplyr::recode(scenario, 
                           "suitability_rcp26" = "rcp26", 
                           "suitability_rcp45" = "rcp45", 
                           "suitability_rcp85" = "rcp85"))

# get individual species and species richness per country: 
df_richness_country_current <- df_compare_future_species %>% 
  dplyr::filter(current == 1) %>% 
  dplyr::select(acc, country_pooled) %>%
  distinct() %>% 
  group_by(country_pooled) %>% 
  mutate(richness_current = n_distinct(acc)) #%>% 
  #mutate(scenario == 'current')

df_richness_country_rcp45 <- 
  df_compare_future_species %>% 
  dplyr::filter(rcp45 == 1) %>% 
  dplyr::select(acc, country_pooled) %>%
  distinct() %>% 
  group_by(country_pooled) %>% 
  mutate(richness_rcp45 = n_distinct(acc)) #%>% 
 # mutate(scenario = 'rcp45')

df_suitab_countries <- df_richness_country_current %>% 
  full_join(df_richness_country_rcp45)


# Calculate frequency of occurrence per species, suitability, and scenario
df_suitability_freq <- df_long_suitability %>%
  group_by(acc, suitability, scenario) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(acc, scenario, desc(count)) %>% 
  dplyr::filter(suitability != 0) %>% 
  mutate(acc = factor(acc),
         suitability = factor(suitability),
         scenario    = factor(scenario   ))


### Sankey plot test: need to restructure data -----------------

df_current <- df_compare_future_species %>% 
  dplyr::filter(current == 1) %>% 
  dplyr::select(site, acc) %>% 
  mutate(scenario = 'current')

df_rcp26 <- df_compare_future_species %>% 
  dplyr::filter(current == 1 & rcp26 == 1) %>% 
  dplyr::select(site, acc) %>% 
  mutate(scenario = 'rcp26')

df_rcp45 <- df_compare_future_species %>% 
  dplyr::filter(current == 1& rcp45 == 1) %>% 
  dplyr::select(site, acc) %>% 
  mutate(scenario = 'rcp45')

df_rcp85 <- df_compare_future_species %>% 
  dplyr::filter(current == 1& rcp85 == 1) %>% 
  dplyr::select(site, acc) %>% 
  mutate(scenario = 'rcp85')

df_rbind = rbind(df_current, df_rcp26, df_rcp45, df_rcp85)

# each species needs to exist across all scenarios
df_alluvial_complete <- df_rbind %>%
  dplyr::filter(acc %in% top_species_site_share) %>% 
  count(site, scenario, acc) %>%  # Count occurrences of each species per site per scenario
  complete(site, scenario, acc, fill = list(n = 0)) %>%  # Ensure all combinations exist
  rename(freq = n) %>%  # Rename count column for clarity
  arrange(site, acc, scenario) #%>% 
#dplyr::filter(freq !=0)

# calculate frequency across species and scenarios
df_alluvial_complete_freq <- df_alluvial_complete %>%
  mutate(scenario = factor(scenario, levels = c("current", "rcp26", "rcp45", "rcp85")),
         acc = factor(acc),
         site = factor(site)) %>% 
  group_by(scenario, acc) %>% 
  dplyr::summarise(sum_n = sum(freq, na.rm = T))

# Define species order based on species_colors
species_order <- names(species_colors)

# Convert acc to a factor with the specified order
df_alluvial_complete_freq$acc <- factor(df_alluvial_complete_freq$acc, levels = species_order)

# Create the alluvial plot with ordered species and custom colors
p.species.clim.suitability <-  ggplot(df_alluvial_complete_freq,
                                      aes(x = scenario, 
                                          y = sum_n,
                                          stratum = acc,
                                          alluvium = acc,
                                          fill = acc)) +  
  geom_flow() +               
  geom_stratum() +  
  theme_classic2(base_size = 8) +
  labs(title = "",
       x = "Scenario",
       y = "Counts of plots with respective species [#]",
       fill = "Species") +
  scale_fill_manual(values = species_colors,
                    labels = species_labels) +
  #scale_fill_manual() +# Replace y-axis labels with full Latin names
  theme(legend.position = "right", #c(0.75,0.7),
    legend.text = element_text(face = "italic", size = 8), # Italicize only legend items
    axis.text.y = element_text(size = 8),  # Keep y-axis labels non-italic
    axis.text.x = element_text(size = 8),  # Adjust x-axis label size
    axis.title = element_text(size = 8)    # Adjust axis title size
  )

p.species.clim.suitability

p.clim.suitab.fin <- ggarrange(p.species.clim.suitability, p.species.suitability_country,nrow = 2, ncol = 1, labels = c("[a]", "[b]") )

# Save the combined plot as an image
ggsave(
  filename = "outFigs/p.clim.suitab.fin.png",    # File name (change extension for different formats, e.g., .pdf)
  plot = p.clim.suitab.fin,             # Plot object to save
  width = 5,                       # Width of the saved plot in inches
  height = 7,                       # Height of the saved plot in inches
  dpi = 300,                        # Resolution (dots per inch)
  units = "in"                      # Units for width and height
)
# p.species.suitability_country_simpl





### Summary table per country -----------------------------


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
plots_no_species <- 
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



# Compare presence with future scenarios: count the number of plots that have not a single species spared between curremt and future scenarios
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

### Barplot of species suitability -------------
# seems wrong ---
df_species_presence <- df_compare_future_species %>%
  pivot_longer(cols = c(current, rcp26, rcp45, rcp85), 
               names_to = "scenario", 
               values_to = "presence") %>%
  filter(presence == 1) %>%  # Keep only species that are present
  select(country_pooled, acc, scenario) %>%
  distinct() %>%  # Remove duplicates
  mutate(presence = acc) %>%  # Rename `acc` column for clarity
  pivot_wider(names_from = scenario, values_from = presence, values_fill = "0") 

# make sankey plot:
# how often which species will transit into each climate change scenarios?



  
  
  
  # Calculate frequency of occurrence grouped by species, scenario, and presence
  df_freq <- df_compare_future_species %>%
    dplyr::filter(current == 1) %>%
    pivot_longer(cols = c(current, rcp26, rcp45, rcp85), names_to = 'scenario', values_to = 'presence') %>%
    mutate(scenario = factor(scenario, levels = c("current", "rcp26", "rcp45", "rcp85")),
           presence = factor(presence)) %>%
    group_by(acc, scenario, presence) %>%
    summarise(freq = n(), .groups = "drop")
  
  

# calculate species richness per country and scanerio - cross link, get a barplot with categories of occurance
# Count unique species in 'current' per country
species_current <- df_compare_future_species %>%
  dplyr::filter(current == 1) %>%  # Filter where species is present
  group_by(country_pooled) %>%
  dplyr::reframe(n_species_current = unique(acc))

# Count unique species in each future scenario per country
#species_scenarios <- 
  df_compare_future_species %>%
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

