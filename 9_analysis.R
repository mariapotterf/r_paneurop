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

# add manually created regions
regions_manual <- terra::vect('rawData/regions_manual2.gpkg')
regions_manual_df <- as.data.frame(regions_manual)

df_fin <- df_fin %>% 
  full_join(regions_manual_df) %>% 
  rename(region_manual = name)

# ad disturbance severity based on presence/absence of mature trees 
df_mature_dist_severity <- fread('outData/disturb_severity_mature.csv')

# Create the opposite vector: listing residual trees
df_mature_dist_severity$residual_mature_trees <- 1 - df_mature_dist_severity$mature_dist_severity


df_fin <- df_fin %>% 
  left_join(df_mature_dist_severity, by = c('site' = 'cluster')) %>% 
  mutate(time_since_disturbance = 2023 - disturbance_year)

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

prop.table(table(df_fin$disturbance_year))

## check for Czechia -------------
df_fin_cz <- df_fin %>% 
  dplyr::filter(country == "CZ")


dim(df_fin_cz)
table(df_fin_cz$dominant_species)

# Assuming the table 'table(df_fin$dominant_species)' has already been created
# as "dominant_species_table"
dom_species_table_cz <- table(df_fin_cz$dominant_species)

# Convert the table to a data frame for easier manipulation
dom_species_table_cz_df <- as.data.frame(dom_species_table_cz)
colnames(dom_species_table_cz_df) <- c("species", "count")

# Calculate proportions
dom_species_table_cz_df$proportion <- round(dom_species_table_cz_df$count / sum(dom_species_table_cz_df$count) * 100, 1)

# Sort the data frame by proportion in descending order
dom_species_table_cz_df <- dom_species_table_cz_df[order(-dom_species_table_cz_df$proportion), ]

# Calculate species richness (number of unique species)
species_richness_cz <- nrow(dom_species_table_cz_df)

# Print the results
print(dom_species_table_cz_df)
cat("Species Richness:", species_richness_cz, "\n")



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


# Assign colors to each species in the order of `top_species_site_share$Species`
species_colors <- setNames(
  reversed_colors,
  c("piab", "fasy", "quro", "pisy", "soau", "acps", "potr", "abal", "besp", "lade")
)

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




## Get summary tables ---------------------------------------------------------

# how many sites do not have any stem density present? not evemn mature treees?
# Count rows and calculate proportions
total_sites <- length(unique(df_fin$site))


# How many plots per dominant tree species ? (rIVI)
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


## get plots share by most frequent dominant species: rIVI -------------------------------
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
    y = "Proportion (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis text since it's a single bar
        axis.title.x = element_blank(),
        
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8, face = "plain")
  )


## share of plots by species richness --------------------------------------------


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

# Create a stacked bar plot
p_bar_richness_groups <- ggplot(richness_summary, aes(x = 1, y = proportion, fill = richness_category)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste0(round(proportion, 1), "%", '(', site_count, ')')), 
            position = position_stack(vjust = 0.5), size = 2) +  # Add percentage labels
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    title = "",
    x = NULL,
    y = "Proportion (%)",
    fill = "Count species\nrichness"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text since it's a single bar
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, face = "plain"),
    plot.title = element_text(size = 8, face = "plain", hjust = 0.5)
  )


## vertical classes share -----------------------------------

# Categorize richness into bins and calculate proportions
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
  geom_bar(stat = "identity", color = "black") +
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
  summarise(sum_stem_density = sum(stem_density, na.rm = t)) #%>% 
  #dplyr::filter(sum_stem_density>0)



## Share of plots by species -----------------------------------------------------

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


### top species CZ start !!!-----

df_stem_sp_sum_cz <- df_stem_sp_sum[grep("15_|26_", df_stem_sp_sum$cluster), ]
# Count the unique clusters (sites) for each species
species_site_share_cz <- df_stem_sp_sum_cz %>%
  dplyr::filter(sum_stem_density >0) %>% 
  group_by(Species) %>%
  summarise(site_count = n_distinct(cluster)) %>%
  ungroup()


# Add a column for the share (percentage) of sites
species_site_share_cz <- species_site_share_cz %>%
  mutate(share_of_sites = (site_count / 135) * 100)

# Sort the results in descending order of share
species_site_share_cz <- species_site_share_cz %>%
  arrange(desc(share_of_sites))

# Print the result
print(species_site_share_cz)

top_species_site_share_cz <- species_site_share_cz%>% 
  arrange(desc(share_of_sites)) %>%
  slice_head(n = 10) # Select the top 7 rows


# !!! cz top end



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
    x = "Plots share\n[%]",
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


### Get a color scheme per species -------------------------

# Reverse the color palette and map to the species in the desired order
n_colors <- 10  # Number of species
my_colors <- colorRampPalette(brewer.pal(11, "RdYlGn"))(n_colors)  # Generate colors

# Reverse the color order to start with dark green
reversed_colors <- rev(my_colors)


# Print the color assignments for confirmation
print(species_colors)


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

#my_species_levels <-  levels(df_stem_sp_sum_ordered$Species)

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
    color = "black",  # Black outline for points
    shape = 21,  # Shape 21 is a circle with a fill and border
    size = 0.5,
    linewidth = 0.2,
    stroke = 0.2,
    position = position_nudge(y = 0.2)  # No vertical nudge for alignment
  ) +
  theme_classic() +
  labs(
    title = "",
    x = "Stem density\n(log10) [#/ha]",
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




### !!!get the species density plot for Czechia -------------

df_cz <- df_stem_sp_sum[grep("15_|26_", df_stem_sp_sum$cluster), ]


# Calculate median for each species and reorder the factor levels
df_stem_sp_sum_ordered <- df_cz %>%
  #dplyr::filter(cou)
  dplyr::filter(sum_stem_density >0) %>% 
  dplyr::filter(Species %in% top_species_site_share_cz$Species ) %>%  #top_species_overall
  dplyr::group_by(Species) %>%
  dplyr::mutate(median_stem_density = median(sum_stem_density, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>%
  mutate(Species = factor(Species, levels = rev(top_species_site_share_cz$Species))) # Set custom order
# dplyr::mutate(Species = reorder(Species, median_stem_density))  # Reorder species by median stem density

# Add a log-transformed column for sum_stem_density
df_stem_sp_sum_ordered <- df_stem_sp_sum_ordered %>%
  mutate(log_sum_stem_density = log10(sum_stem_density + 1))  # Adding 1 to avoid log(0)


df_stem_sp_sum_ordered <- df_stem_sp_sum_ordered %>% 
  group_by(Species) %>% 
  mutate(median_value = median(log_sum_stem_density, na.rm = TRUE)) %>% 
  arrange(desc(median_value)) %>% 
  ungroup()



p_stem_density_species_cz <- df_stem_sp_sum_ordered %>%
  ggplot(aes(x = log_sum_stem_density, y = Species, group = Species)) +
  geom_density_ridges(
    aes(fill = Species), 
    alpha = 1, 
    color = 'NA', 
    scale = 0.9 # Adjust the vertical scale of the ridges
  ) +
 # scale_fill_manual(values = species_colors) +
  stat_summary(
    aes(x = log_sum_stem_density, fill = Species),  # Add fill aesthetic for inner color
    fun = median, 
    fun.min = function(x) quantile(x, 0.25),  # 25th percentile (Q1)
    fun.max = function(x) quantile(x, 0.75),  # 75th percentile (Q3)
    geom = "pointrange", 
    color = "black",  # Black outline for points
    shape = 21,  # Shape 21 is a circle with a fill and border
    size = 0.5,
    linewidth = 0.2,
    stroke = 0.2,
    position = position_nudge(y = 0.2)  # No vertical nudge for alignment
  ) +
  theme_classic() +
  labs(
    title = "",
    x = "Stem density\n(log10) [#/ha]",
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
print(p_stem_density_species_cz)

# CZEchia stem density end ----



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
  summarise(med_share = median(share))

# Transform data to make saplings negative for diverging plot
species_composition_plot_data <- species_composition_sapl_juv_med %>%
  mutate(
    med_share = ifelse(VegType == "Saplings", -med_share, med_share) # Saplings are negative
  ) %>% 
  mutate(Species = factor(Species, 
                          levels = rev(top_species_site_share$Species))) # Set custom order

  
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

# Diverging bar chart with IQR
ggplot(species_composition_sapl_juv_summary, aes(x = med_share, y = Species, fill = VegType)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.8) + # Horizontal bars
  geom_errorbarh(aes(xmin = iqr_low, xmax = iqr_high), height = 0.2, color = "black") + # Add IQR error bars
  scale_x_continuous(
    labels = abs, # Show positive labels
    name = "Median Share (%)"
  ) +
  scale_fill_manual(values = species_colors) + # Apply custom color palette
  labs(
    x = "",
    y = "Species",
    title = "Diverging Bar Chart with IQR: Share of Saplings and Juveniles"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, hjust = 0.5)
  )


## adjust colors -----------------------------------------------------------------
library(scales) # For color adjustment (lighter/darker shades)

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
p_share_vertical_species <- ggplot(species_composition_sapl_juv_summary, aes(x = med_share, y = Species, fill = fill_key)) +
  geom_bar(stat = "identity", position = "identity", alpha = 1) + # Horizontal bars
  geom_errorbarh(aes(xmin = iqr_low, xmax = iqr_high), 
                 height = 0.1, 
                 color = "black",
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
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.2) + #
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


library(cowplot)

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

# Remove x- and y-axis labels from individual plots
p1 <- p_sites_share_species +
  scale_y_discrete(labels = species_labels) + # Replace y-axis labels with full Latin names
 # theme_classic(base_size = 8) +
  theme(
    axis.text.y = element_text(face = "italic", size = 8) #, # Italics for Latin names
 #   axis.text.x = element_text(size = 8),  # Adjust x-axis label size
  #  axis.title = element_text(size = 8)   # Adjust axis title size
  )
p2 <- p_stem_density_species + theme(axis.ticks.y = element_blank(),
                                     axis.text.y = element_blank(), axis.title.y = element_blank())
p3 <- p_share_vertical_species + theme(axis.ticks.y = element_blank(), 
                                       axis.text.y = element_blank(), axis.title.y = element_blank())


# Combine plots with aligned axes
combined_plot <- plot_grid(
  p1,p2,p3,
  align = "h",  # Align vertically
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
  dplyr::filter(Species %in% top_species_site_share$Species) %>% 
  dplyr::filter(stem_density > 0) %>% 
  mutate(Species = factor(Species, levels = rev(my_species_levels ))) %>%  # Set ordertop_species_site_share$Species
  mutate(
   # Species = factor(Species, levels = top_species_site_share$Species),  # Set order
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
  dplyr::filter(Species %in% top_species_site_share$Species) %>% 
  dplyr::filter(stem_density > 0) %>% 
  mutate(Species = factor(Species, levels = my_species_levels)) %>%  # Set order
  mutate(
    Species = factor(Species, levels = top_species_site_share$Species),  # Set order
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






## Tableplot:  Stem density per species and vertical class: --------------------------------

# Get a summary table:
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
  dplyr::filter(Species %in% top_species_site_share$Species  ) %>% 
  dplyr::filter(stem_density > 0) %>% 
  mutate(Species = factor(Species, 
                          levels = top_species_site_share$Species)) %>%
  
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
species_composition <- stem_dens_species_long_cluster %>%
  left_join(df_fin_country_ind, by = c("cluster" = "site")) %>% 
  group_by(Species, country_abbr) %>%
  summarize(sum_stems = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the total stem density per climate class and the share of each species
species_composition <- species_composition %>%
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
p_species_distribution_country <- species_composition %>% 
  dplyr::filter(share >5) %>% 
  ggplot(
                                 aes(x = country_abbr, 
                                     y = share, 
                                     fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  geom_text(aes(label = ifelse(share >= 5, paste0(round(share, 1), ""), "")),
            position = position_stack(vjust = 0.5),  # Labels inside the bars
            size = 3, color = "black") +  # Adjust text size and color
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
df_fin$sum_stems_mature_pres_abs <- ifelse(df_fin$sum_stems_mature > 0, 1, 0)





#### check for multicollinearity -----------------------------------------------------

library(car)
selected_data <- df_fin %>%
  dplyr::select(stem_regeneration, 
                drought_spei12,
                spei12,
                tmp,  
                prcp,
                spei1,
                cv_t2m,
                cv_tp,
                sd_grw_anm_tmp
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
m1 <- gam(stem_density ~ s(spei12, k = 15),  # Factors included without s() for categorical variables
          family = nb,  # Negative Binomial to handle overdispersion
          data = df_fin)

# tw distribution is the best


# TW has a better fit, also can handle zero!
# select main variables as predictors 
predictor_vars_sub <- c(
  "spei12",
  "tmp", 
  "prcp", 
  
  "drought_tmp",
  "drought_prcp",
  
  "tmp_z", 
  "prcp_z", 
  
  # get drorught indiced : spei in 2018-2019
  "drought_spei1",
  "drought_spei12",
  
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
   #"sum_stems_mature_scaled" , # mature trees stems: scaled
   "residual_mature_trees",        
   #"sum_stems_mature_pres_abs",  # mature trees present/absent
  
  "clay_extract", 
  "depth_extract", 
   
   # site info
  "av.nitro",
  
  # anomalies per growing season for tmp and prcp
  "mean_grw_anm_prcp" ,        
  "mean_grw_anm_tmp"  ,        
  "sd_grw_anm_prcp",           
  "sd_grw_anm_tmp" ,           
  "median_grw_anm_prcp"   ,   
  "median_grw_anm_tmp" ,      
  "max_grw_anm_prcp" ,        
  "max_grw_anm_tmp",
 # "min_grw_anm_prcp",
#  "min_grw_anm_tmp",  
  
  # seasonality: CV - over year
  "cv_t2m",
  "cv_tp"    
  #"richness",
  #'rIVI',
  
  #  "n_vertical"
  )



# run univariate models to find teh best predictors :

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
#pairs(eda_data, main = "")



## Models: simplify table just for stem_regeneration  ---------------------------------------------------------------------------------
# test drivers: simplify the analysis:
# Subset the data

df_stem_regeneration2 <- df_fin %>% 
  dplyr::select(all_of(c("site", "stem_regeneration", predictor_vars_sub,
                         
                         "country_pooled","region", "region_manual", "clim_grid",  "x", "y")))


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
 # dplyr::select(where(is.numeric))
  dplyr::select(prcp, tmp, spei12,
                drought_prcp, drought_tmp, drought_spei12,
                sd_grw_anm_prcp,
                sd_grw_anm_tmp,
                tmp_z, prcp_z,
                cv_t2m,
                cv_tp,
                distance_edge, depth_extract,
                disturbance_severity, #mature_dist_severity ,
                #sum_stems_mature ,
                depth_extract , clay_extract, av.nitro, #management_intensity
                )

# Calculate the correlation matrix
correlation_matrix <- cor(predictors, use = "complete.obs")

# Display the correlation matrix
library(corrplot)

# Plot using corrplot
corrplot(correlation_matrix, method = "color", type = "upper",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         #addCoef.col = "black", 
         tl.col = "black", tl.srt = 45)

par(mfrow = c(1, 2))
plot(df_fin$prcp, df_fin$clim_grid)
plot(df_fin$tmp, df_fin$clim_grid )
dev.off()

# understand teh one to one relationship --------------

ggplot(df_stem_regeneration2, aes(x = drought_spei12, y = stem_regeneration     )) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue")

ggplot(df_stem_regeneration2, aes(x = sd_grw_anm_tmp, y = stem_regeneration     )) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue")

ggplot(df_stem_regeneration2, aes(x = clay_extract, y = stem_regeneration     )) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue")

ggplot(df_stem_regeneration2, aes(x = cv_t2m, y = stem_regeneration     )) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue")

# tmp_z is anomaly over whole year
ggplot(df_stem_regeneration2, aes(x = tmp_z, y = stem_regeneration     )) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue")


ggplot(df_stem_regeneration2, aes(x = max_grw_anm_tmp, y = stem_regeneration     )) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue")

ggplot(df_stem_regeneration2, aes(x = median_grw_anm_prcp, y = stem_regeneration     )) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue")

# get stats: distnace to edge and severity typesper country
ggplot(df_stem_regeneration2, aes(x = country_pooled , y = distance_edge     )) +
  geom_boxplot()
  
# get stats: distnace to edge and severity typesper country
ggplot(df_stem_regeneration2, aes(x = country_pooled , y = disturbance_severity     )) +
  geom_boxplot()

# get stats: distnace to edge and severity typesper country
ggplot(df_stem_regeneration2, aes(x = country_pooled , y = mature_dist_severity     )) +
  geom_boxplot()



# Export as a PNG
png("outFigs/correlation_matrix_plot.png", width = 800, height = 800, res = 300)
plot.new()
corrplot(correlation_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, 
         title = "Correlation Matrix of Predictors",
         addCoef.col = "black", number.cex = 0.7)
dev.off()#summary(interaction_model_4)




# et median data per clim cluster -----------------------------------------one for climate: on grid 9x9, for disturbance severity vs distance_edge - on local basis

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
               #s(spei12, k = 5) + 
               #s(distance_edge, k = 5) +
               #s(depth_extract, k = 4) +
               #s(disturbance_severity, k =5) +
               #s(mature_dist_severity, k = 5) +
               #s(clay_extract, k = 5) +
               #s(av.nitro, k =5) +
               ti(tmp, prcp, k = 5) +
               te(disturbance_severity, distance_edge, k = 5) +
               #ti(disturbance_severity_c, prcp_c, k = 5) +
               #s(management_intensity,by = country_pooled, k = 4) + 
               #s(country_pooled, bs = "re") + 
               s(x,y) #+
               #s(clim_grid, bs = "re") 
             ,
             family = tw(), method = "REML", data = df_median)


### the most meaningful: use disturbance severity, distance to edge ----------------------
m_int_sev_edge <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    ti(disturbance_severity,distance_edge, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    #s(country_pooled, bs = "re") +
    #s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  data = df_stem_regeneration2
)

cor(df_stem_regeneration2$disturbance_severity, df_stem_regeneration2$residual_mature_trees) # high correlation, keep only one of them!



# only interaction modelling: for climate and disturbance chars 
m_int_res_edge_full_te <- gam(
  stem_regeneration ~ 
    #s(prcp, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(residual_mature_trees, distance_edge, k = 5 ) +
    te(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    #s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


# interaction and single effects for climate
# just interaction modelling
m_int_res_edge_full_te_comb <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(residual_mature_trees, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    #s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)



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

# add seasonality
m_int_sev_edge_full_te_comb2 <- gam(
  stem_regeneration ~  s(drought_prcp, k = 5) +
    s(drought_tmp, k = 5) +
    s(sd_grw_anm_tmp , k = 5) +
    s(av.nitro, k = 5) +
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


m_int_sev_edge_full_te_comb3 <- gam(
  stem_regeneration ~  s(drought_prcp, k = 5) +
    s(drought_tmp, k = 5) +
    s(sd_grw_anm_tmp , k = 5) +
    s(av.nitro, k = 5) +
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp, tmp, k = 5 ) +
    ti(drought_prcp,sd_grw_anm_tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)
p<- ggpredict(m_int_sev_edge_full_te_comb3, )

AIC(m_int_sev_edge_full_te_comb2, m_int_sev_edge_full_te_comb, m_int_sev_edge_full_te_comb3)
summary(m_int_sev_edge_full_te_comb3)

m_int_spei1 <- gam(
  stem_regeneration ~ 
    s(drought_spei1, k = 5) + s(tmp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(drought_spei1,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m_int_spei12 <- gam(
  stem_regeneration ~ 
    s(drought_spei12, k = 5) + s(tmp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(drought_spei12,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)



m_int_sev_edge_full_te_comb_factors <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    country_pooled +
    #region_manual +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)
AIC(m_int_sev_edge_full_te_comb_factors, m_int_sev_edge_full_te_comb)
summary(m_int_sev_edge_full_te_comb)
summary(m_int_sev_edge_full_te_comb_factors)


### 20250228 add seasonality: CV/sd ------------------------------------------------

m_sd_tmp <- gam(
  stem_regeneration ~  s(sd_grw_anm_tmp, k = 5) +
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

m_sd <- gam(
  stem_regeneration ~  s(sd_grw_anm_tmp, k = 5) +
    s(prcp, k = 5) +# s(tmp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,sd_grw_anm_tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)
summary(m_sd)

# start simple with CVS
m1 <- gam(
  stem_regeneration ~  s(sd_grw_anm_tmp, k = 5) +
    s(drought_spei1, k = 5) +
    s(prcp, k = 5) + 
    s(tmp, k = 5) +
    s(cv_t2m, k = 5) + 
    s(cv_tp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
   # te(disturbance_severity, distance_edge, k = 5 ) +
    #ti(prcp,sd_grw_anm_tmp, k = 5 ) +
    #s(management_intensity,by = country_pooled, k = 4) + 
    #s(country_pooled, bs = "re") +
    #s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


m2 <- gam(
  stem_regeneration ~  s(sd_grw_anm_tmp, k = 5) +
    s(drought_spei1, k = 5) +
    s(prcp, k = 5) + 
    s(tmp, k = 5) +
    s(cv_t2m, k = 5) + 
    s(cv_tp, k = 5) +
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

m3 <- gam(
  stem_regeneration ~  s(sd_grw_anm_tmp, k = 5) +
    s(drought_spei12, k = 5) +
    s(prcp, k = 5) + 
    s(tmp, k = 5) +
    s(av.nitro, k = 5) +
    #s(cv_t2m, k = 5) + 
    s(cv_tp, k = 5) +
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


m4 <- gam(
  stem_regeneration ~  s(sd_grw_anm_tmp, k = 5) +
    s(drought_prcp, k = 5) +
    s(prcp, k = 5) + 
    s(tmp, k = 5) +
    s(av.nitro, k = 5) +
    #s(cv_t2m, k = 5) + 
    s(cv_tp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
   # ti(drought_prcp,sd_grw_anm_tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


m5 <- gam(
  stem_regeneration ~  s(sd_grw_anm_tmp, k = 5) +
    s(drought_prcp, k = 5) +
    #s(prcp, k = 5) + 
   # s(tmp, k = 5) +
    s(av.nitro, k = 5) +
    #s(cv_t2m, k = 5) + 
    s(cv_tp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    #te(prcp,tmp, k = 5 ) +
    ti(drought_prcp,sd_grw_anm_tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)
summary(m5)

# Interaction 1: Precipitation and Temperature
p1 <- ggpredict(m_int_sev_edge_full_te_comb3, terms = c("drought_prcp"))
plot(p1)

p1 <- ggpredict(m_int_sev_edge_full_te_comb3, terms = c("sd_grw_anm_tmp"))
plot(p1)


p1 <- ggpredict(m_int_sev_edge_full_te_comb3, terms = c("drought_prcp", "sd_grw_anm_tmp[1,1.4]"))
plot(p1)


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
             size = 1, alpha = 0.2, color = 'grey') +
  geom_line(linewidth = 1, aes(color = group) ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000, fill = group), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  scale_fill_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  theme_classic() +
  labs(x = "Precipitation [mm]", 
       y = "Regeneration stem density [#*1000/ha]", title = "p<0.0001", 
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






# check VIF again
selected_data <- df_fin %>% dplyr::select(stem_regeneration,  prcp,  sd_grw_anm_tmp, cv_t2m, cv_tp)# drought_spei1,tmp,
lm_model <- lm(stem_regeneration ~ ., data = selected_data)
vif(lm_model)
# still low VIF

final_model <- gam(stem_density ~ s(drought_spei1, k = 5) + 
                     s(prcp, k = 5) + 
                     s(tmp, k = 5) + 
                     s(sd_grw_anm_tmp, k = 5) + 
                     s(cv_t2m, k = 5) + 
                     s(cv_tp, k = 5),
                   family = tw(), data = df_fin)

AIC(final_model, m_sd_tmp, m_int_sev_edge_full_te_comb, m2, m3)

# store the best model for regeneration density
#fin.m.reg.density <- m_int_sev_edge_full_te_comb     

fin.m.reg.density <- m_int_sev_edge_full_te_comb#m_anom_prcp_re_reg

vis.gam(fin.m.reg.density, view = c("prcp", "tmp"), plot.type = "persp",
        main = "Interaction between Precipitation and Temperature",
        zlab = "Stem Regeneration", xlab = "Precipitation", ylab = "Temperature")


cor_matrix <- cor(df_fin %>% dplyr::select(drought_spei1, drought_spei12, tmp, prcp), use = "pairwise.complete.obs", method = "pearson")
(cor_matrix)
 
# drought_spei1 drought_spei12        tmp       prcp
# drought_spei1      1.0000000      0.7803899 -0.1644824  0.6986090
# drought_spei12     0.7803899      1.0000000 -0.1696574  0.7655749
# tmp               -0.1644824     -0.1696574  1.0000000 -0.4270087
# prcp               0.6986090      0.7655749 -0.4270087  1.0000000


### 20250220 - ry to improve teh best mdoel with the updated DWD values ---------
m_int_sev_edge_full_te_comb_spei <- gam(
  stem_regeneration ~ s(spei12, k = 5) + 
    #s(prcp, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    # s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m_int_sev_edge_full_te_comb_spei_add <- gam(
  stem_regeneration ~ s(spei12, k = 5) + 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    # s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


m_int_sev_edge_full_te_comb_spei_add2 <- gam(
  stem_regeneration ~ s(spei12, k = 5) + 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(spei12, prcp, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    # s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


m_int_sev_edge_full_te_comb_spei_add2_re <- gam(
  stem_regeneration ~ s(spei12, k = 5) + 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(spei12, prcp, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


m_int_sev_edge_full_te_comb_spei_add2_re_simpl <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    #ti(spei12, prcp, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

### anomalies -------------------
m_anom_tmp_z <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(prcp, k = 5) + s(tmp_z, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    #ti(spei12, prcp, k = 5 ) +
    ti(prcp,tmp_z, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
   # s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m_anom_prcp_z <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(prcp_z, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    #ti(spei12, prcp, k = 5 ) +
    ti(prcp_z,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    # s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m_anom_prcp_z_re <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(prcp_z, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    #ti(spei12, prcp, k = 5 ) +
    ti(prcp_z,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(clim_grid, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


m_anom_prcp_z_re_reg <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(prcp_z, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    #ti(spei12, prcp, k = 5 ) +
    ti(prcp_z,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region, bs = "re", k = 5) +                # Macro-scale random effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


m_anom_prcp_re_reg <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    #ti(spei12, prcp, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # simply the macro scale effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m_anom_prcp_fix_reg <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(residual_mature_trees  , k = 5) + 
    #s(sum_stems_mature  , k = 5) + 
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    #ti(spei12, prcp, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    region_manual +                # use as fixed effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

# tyry again to include teh dist values as single factors along intercation:

m_anom_prcp_re_reg_test <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = "re", k = 5) +                # simply the macro scale effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m_anom_prcp_re_reg_test_fix <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    region_manual +                # simply the macro scale effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


AIC(m_anom_prcp_re_reg_test, m_anom_prcp_re_reg, m_anom_prcp_re_reg_test_fix)



AIC(m_int_sev_edge_full_te_comb,
    m_int_sev_edge_full_te_comb_spei, 
    m_int_sev_edge_full_te_comb_spei_add2, 
    m_int_sev_edge_full_te_comb_spei_add2_re,
    m_int_sev_edge_full_te_comb_spei_add2_re_simpl,
    m_anom_tmp_z,
    m_anom_prcp_z,
    m_anom_prcp_z_re,
    m_anom_prcp_z_re_reg,
    m_anom_prcp_fix_reg)



### 20250227 check anomalies over growth season --------------------------------
# more difficult to interpret - skip 
m0 <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(tmp, k = 5) + s(prcp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(tmp,prcp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = 're') +                # simply the macro scale effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


m1 <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(max_grw_anm_tmp, k = 5) + s(max_grw_anm_prcp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(max_grw_anm_prcp,max_grw_anm_tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = 're') +                # simply the macro scale effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m1_te <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(max_grw_anm_tmp, k = 5) + s(max_grw_anm_prcp, k = 5) +
   # s(distance_edge, k = 5) +
   # s(disturbance_severity, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(max_grw_anm_prcp,max_grw_anm_tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = 're') +                # simply the macro scale effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m2 <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(tmp_z, k = 5) + s(prcp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(tmp_z,prcp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = 're') +                # simply the macro scale effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)

m3 <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(sd_grw_anm_tmp, k = 5) + s(max_grw_anm_prcp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(max_grw_anm_prcp,sd_grw_anm_tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = 're') +                # simply the macro scale effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)


m4 <- gam(
  stem_regeneration ~ #s(spei12, k = 5) + 
    s(median_grw_anm_tmp, k = 5) + s(median_grw_anm_prcp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(median_grw_anm_prcp,median_grw_anm_tmp, k = 5 ) +
    s(management_intensity,by = country_pooled, k = 4) + 
    s(country_pooled, bs = "re") +
    s(region_manual, bs = 're') +                # simply the macro scale effect
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration2
)
AIC(m0, m1,m1_te, m2,m3, m4)
summary(m1_te)
summary(m4)

fin.m <- m1

# how does temp and max_growth tmp anomalies looks like ?

df_stem_regeneration2 %>% 
  ggplot(aes(x = median_grw_anm_tmp,
             y = tmp)) + 
  geom_point()

 
# plot in easy way how tdoes the k value affect interaction 
p1 <- ggpredict(fin.m, terms = "max_grw_anm_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m, terms = "max_grw_anm_prcp [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m, terms = c("max_grw_anm_prcp", "max_grw_anm_tmp[3,4]"), allow.new.levels = TRUE)
p4 <- ggpredict(fin.m, terms = c("distance_edge", "disturbance_severity[0.3,0.9]"), allow.new.levels = TRUE)



# test simple plot:
ggplot(p3, aes(x = x , y = predicted , ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) 





### check variability withing clim grid --------------------------

# Calculate standard deviation of stem_regeneration for each clim_grid
regional_variability <- df_fin %>%
  dplyr::filter(region_manual != 'CH2') %>% 
  group_by(region_manual) %>%
  summarise(
    mean_stem_regeneration = mean(stem_regeneration, na.rm = TRUE),
    sd_stem_regeneration = sd(stem_regeneration, na.rm = TRUE),
    cv_stem_regeneration = (sd_stem_regeneration / mean_stem_regeneration) * 100,
    
    mean_elevation = mean(elevation, na.rm = TRUE),
    sd_elevation = sd(elevation, na.rm = TRUE),
    cv_elevation = (sd_elevation / mean_elevation) * 100,
    
    mean_tmp = mean(tmp, na.rm = TRUE),
    sd_tmp = sd(tmp, na.rm = TRUE),
    cv_tmp = (sd_tmp / mean_tmp) * 100,
    
    mean_prcp = mean(prcp, na.rm = TRUE),
    sd_prcp = sd(prcp, na.rm = TRUE),
    cv_prcp = (sd_prcp / mean_prcp) * 100
  )

# View variability
#View(regional_variability)

# Function to plot coefficient of variation
plot_cv <- function(data, variable, title, y_label) {
  ggplot(data, aes(x = region_manual, y = !!sym(variable))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = title, x = "Region", y = y_label) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Generate plots
p_stems_cv <- plot_cv(regional_variability, "cv_stem_regeneration", "Variability in Stem Regeneration by Region", "CV of Stem Regeneration")
p_tmp_cv <- plot_cv(regional_variability, "cv_tmp", "Variability in Temperature by Region", "CV of Temperature")
p_prcp_cv <- plot_cv(regional_variability, "cv_prcp", "Variability in Precipitation by Region", "CV of Precipitation")
p_elevation_cv <- plot_cv(regional_variability, "cv_elevation", "Variability in Elevation by Region", "CV of Elevation")

ggarrange(p_stems_cv, p_tmp_cv, p_prcp_cv,p_elevation_cv, ncol = 2, nrow = 2)

##### sites per clim_grid ---------------------------------------

site_clim_grid <- df_fin %>%
  group_by(clim_grid, country_abbr) %>%
  summarise(num_sites = n())

summary(site_clim_grid$num_sites)


# get summary statistics of variation in stem density, elevation, tmp and prcp within manual clusters



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


# Plot: Drivers  ---------------------------------------------------------------------------
y_lab = 'Stem density [#/ha]'

summary(fin.m.reg.density)

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


###  dynamic plot title: ad significance level---------------------------------------------

# Function to extract p-values for smooth terms
extract_p_values <- function(model) {
  summary_table <- summary(model)$s.table
  p_values <- summary_table[, "p-value"]
  names(p_values) <- rownames(summary_table) # Assign variable names
  return(p_values)
}

# Function to format p-values into significance levels or as numerical values
format_p_value <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Extract and format p-values from the model
p_values <- extract_p_values(fin.m.reg.density)
formatted_p_values <- sapply(p_values, format_p_value)

# Function to dynamically create plot titles (without s() prefix)
create_dynamic_plot_title <- function(variable, formatted_p_values, model_p_values) {
  # Remove the s() prefix if present
  clean_variable <- gsub("^s\\((.*)\\)$", "\\1", variable)
  
  # Check if the variable exists in formatted_p_values
  if (!is.na(formatted_p_values[variable])) {
    return(paste0(clean_variable, " ", formatted_p_values[variable], 
                   round(model_p_values[variable], 3) ))
  } else {
    return(clean_variable)
  }
}

# Make titles with names

title_disturbance_severity  = create_dynamic_plot_title("s(disturbance_severity)", formatted_p_values, p_values)
title_residual_mature_trees = create_dynamic_plot_title("s(residual_mature_trees)", formatted_p_values, p_values)
title_distance_edge         = create_dynamic_plot_title("s(distance_edge)", formatted_p_values, p_values)
title_interaction1          = create_dynamic_plot_title("ti(prcp,tmp)", formatted_p_values, p_values)




# make plot manually:  --------------------------

m <- fin.m.reg.density #m_int_sev_edge_full_te_comb # m_int_res_edge_full_te_comb #  m_int_res_edge_full_te
k.check(m)
summary(m)


# test -----------
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
              size = 1, alpha = 0.2, color = 'grey') +
  geom_line(linewidth = 1, aes(color = group) ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000, fill = group), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  scale_fill_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  theme_classic() +
  labs(x = "Precipitation [mm]", 
       y = "Regeneration stem density [#*1000/ha]", title = "p<0.0001", 
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
              height = 1) +
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
  labs(x = "Distance to edge [m]", y = "", title = "p=0.001") +
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

summary(df_fin$elevation)


### p combined no points --------------------------

# Plot the first interaction
p1 <- 
  ggplot(pred1, aes(x = x, y = predicted/1000)) +
  # geom_point(data = filtered_df_plot99, 
  #            aes( x = prcp, y = stem_regeneration/1000), 
  #            size = 1, alpha = 0.2, color = 'grey') +
  geom_line(linewidth = 1, aes(color = group) ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000, fill = group), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  scale_fill_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  theme_classic() +
  labs(x = "Precipitation [mm]", 
       y = "Regeneration stem density [#*1000/ha]", title = "p<0.0001", 
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
  # geom_jitter(data = filtered_df_plot90, 
  #             aes( x = distance_edge, y = stem_regeneration/1000), 
  #             size = 1, alpha = 0.2, color = 'grey',
  #             width = 7,
  #             height = 1) +
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
  labs(x = "Distance to edge [m]", y = "", title = "p=0.001") +
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

          
          
# check difference between groups : resiidual mature tree cover    ----------------- 
# presence - absence, to have enought samples in each group
df_fin %>% 
ggplot(aes(x = adv_delayed,
           y = residual_mature_trees)) +
  geom_boxplot()

plot(df_fin$disturbance_severity, df_fin$residual_mature_trees)

table_data <- table(df_fin$adv_delayed, df_fin$sum_stems_mature_pres_abs )


# Display the table in percentages
percent_table <- prop.table(table(df_fin$adv_delayed, df_fin$residual_mature_trees))*100

# Combine levels 0.25 and 0.6 into "Other"
table_data_combined <- table_data
colnames(table_data_combined)[2:5] <- c("Other", "Other", "Other", "Other")
table_data_combined <- rowsum(table_data_combined, group = colnames(table_data_combined))
fisher_test_combined <- fisher.test(table_data_combined)
fisher_test_combined



# export final drivers model -----------------------------------------------------------------

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
  "spei12",
  "drought_spei12",
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
print(summary_stats_narrow, n=20)

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
    geom_errorbar(aes(ymin = Median - IQR / 2, ymax = Median + IQR / 2), 
                  width = 0.2, position = position_dodge(width = 0.5)) +  # Error bars for IQR
    labs(title = "Median and IQR Plot for Delayed vs Advanced Regeneration",
         x = "Variables", y = "Median and IQR") +
    theme_classic() +
    scale_y_continuous(labels = scales::comma) +  # Optional: Make y-axis readable
    facet_wrap(. ~ Variable, scales = 'free') +   # Free scales for different variables
    theme(axis.text.x = element_blank(),          # Hide x-axis text to avoid clutter
          axis.ticks.x = element_blank(),
          legend.title = element_blank())        # Remove legend title

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
                max_grw_anm_prcp,
                max_grw_anm_tmp,
                
                cv_t2m, 
                cv_tp,
                
                #spei12,
                #spei1,
                #spei6,
                drought_spei1,
                drought_spei12,
                drought_tmp,
                drought_prcp,
                disturbance_severity, 
                distance_edge, 
                adv_delayed) %>% 
  #dplyr::select(all_of(variables_to_plot), adv_delayed) %>% 
  gather(key = "Variable", value = "Value", -adv_delayed) %>% 
  droplevels()
unique(df_long_narrow$Variable)

# list groups to pairwise comparison
comparisons <- list(c("Delayed", "Other"), c("Delayed", "Advanced"), c("Other", "Advanced"))


# Compute Wilcoxon effect sizes for each variable
library(rstatix)

effect_sizes <- df_long_narrow %>%
  group_by(Variable) %>%
  wilcox_effsize(Value ~ adv_delayed, comparisons = comparisons)

# most of my effect sizes (differences between two groups) are rather small, the highest effect has drought_prcp
# boot strap analysis: 

library(boot)
boot_effsize <- boot(data = df_long_narrow, statistic = function(data, indices) {
  sample_data <- data[indices, ]
  wilcox_effsize(sample_data$Value ~ sample_data$adv_delayed)
}, R = 1000)


# Plot using ggboxplot
p_boxplot_wilcox_narrow <- ggboxplot(df_long_narrow, x = "adv_delayed", y = "Value", 
                              fill = "adv_delayed", 
                              palette = c("lightblue", "red", "green"),
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
    text = element_text(size = 5),         # Set all text to size 3
    axis.text = element_text(size = 5),    # Axis tick labels
    axis.title = element_text(size = 5),   # Axis titles
    strip.text = element_text(size = 5),   # Facet labels
    legend.text = element_text(size = 5),  # Legend text
    plot.title = element_text(size = 5),   # Plot title
    strip.background = element_blank(),    # Remove the box around facet names
    strip.placement = "outside",           # Optional: Move facet label outside the plot area
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)  # Add a square border around the plots
  )  +
#  p_boxplot_wilcox_narrow + 
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", 
                     size = 2, label.x = 1.5) +
  geom_text(data = effect_sizes, aes(x = 2, y = max(df_long_narrow$Value), 
                                     label = paste0("r = ", round(effsize, 2))))

p_boxplot_wilcox_narrow

# use bayes statistics: more robust for uneven sample sizes
library("BayesFactor")
# Compute Bayesian ANOVA for each variable
bayesian_results <- df_long_narrow %>%
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


# filter out outliers
df_long_narrow_filtered <- df_long_narrow %>%
  group_by(adv_delayed, Variable) %>%
  dplyr::filter(
    Value > quantile(Value, 0.25) - 1.5 * IQR(Value) & 
      Value < quantile(Value, 0.75) + 1.5 * IQR(Value)
  ) %>%
  ungroup()


  #count(adv_delayed)

p_boxplot_wilcox_outliers <- 
 #df_long_narrow_filtered %>% 
  df_long_narrow %>% 
  mutate(Variable = factor(Variable, 
                           levels = c('prcp', 'tmp', "distance_edge", "disturbance_severity",
                                       
                                      #"max_grw_anm_prcp",
                                      #"max_grw_anm_tmp",
                                      "drought_spei1"#,
                                     # "mature_dist_severity", 
                                     # "sum_stems_mature",  
                                     # "residual_mature_trees"
                                     ))) %>% 
  droplevels() %>% 
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
   # geom_jitter(size = 0.1, alpha = 0.5, aes(color = adv_delayed)) +
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
p_boxplot_wilcox_outliers


## START
df_long_narrow_sub <- df_long_narrow %>%
  dplyr::filter(Variable %in% c('prcp', 'drought_spei1', "disturbance_severity")) %>% 
  droplevels() %>% 
  mutate(Variable = factor(Variable, 
                         levels = c("prcp", "disturbance_severity", "drought_spei1"), 
                         labels = c("Precipitation (mm)", "Disturbance Severity", "Drought Index (SPEI)")))

# Define custom labels for facets
facet_labels <- c(
  "prcp" = "Precipitation (mm)", 
  "drought_spei1" = "Drought Index (SPEI)", 
  "disturbance_severity" = "Disturbance Severity"
)

# Create the plot using ggplot instead of ggboxplot
p_boxplot_wilcox_outliers_signif <- df_long_narrow_sub %>% 
  #dplyr::filter(Variable %in% c('prcp', 'drought_spei1', "disturbance_severity")) %>% 
  #droplevels() %>% 
 # mutate(Variable = factor(Variable, levels = c('prcp', 'disturbance_severity', 'drought_spei1'))) %>%
  
  ggplot(aes(x = adv_delayed, y = Value, fill = adv_delayed)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, size = 0.2) +
  scale_fill_manual(values = c("#A50026", "#FDAE61", "#006837")) +
  
  # Add facet labels using facet_wrap
  facet_wrap(~ Variable, scales = "free_y") + #labeller = as_labeller(facet_labels)
  
  # Wilcoxon test annotations
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", size = 2, label.x = 1.5) +
  
  # Titles and axis labels
  labs(title = "", x = "Regeneration Status", y = "Values") +
  
  # Theme modifications
  theme(
    legend.position = 'none',
    text = element_text(size = 5),         
    axis.text = element_text(size = 5),    
    axis.title = element_text(size = 5),   
    strip.text = element_text(size = 5),   # Ensures facet labels show properly
    legend.text = element_text(size = 5),  
    plot.title = element_text(size = 5),   
    strip.background = element_blank(),    
    strip.placement = "outside",           
    panel.grid = element_blank(),          
    axis.line = element_line(color = "black", linewidth = 0.5)  
  )

# Display the plot
print(p_boxplot_wilcox_outliers_signif)
my_colors <- c("#E67E22", "#66A060", "#29502A")
#my_colors <- viridis(3, option = "C")
### END

p_boxplot_wilcox_outliers_signif <-   df_long_narrow %>% 
  dplyr::filter(Variable %in% c('prcp', 'drought_spei1',"disturbance_severity" )) %>% 
  droplevels() %>% 
  mutate(Variable = factor(Variable, 
                           levels = c('prcp', #'tmp', 
                                      #"distance_edge", 
                                      "disturbance_severity",
                                      #"max_grw_anm_prcp",
                                      #"max_grw_anm_tmp",
                                      "drought_spei1"#,
                                      # "mature_dist_severity", 
                                      # "sum_stems_mature",  
                                      # "residual_mature_trees"
                           ))) %>% 
  
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
  # geom_jitter(size = 0.1, alpha = 0.5, aes(color = adv_delayed)) +
  scale_color_manual(values = c("#A50026", 
                                "#FDAE61",
                                "#006837")) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", 
                     size = 2,
                     label.x = 1.5) +  # Position labels between the groups
  labs(title = "",
       x = "", y = "Vals") +
  # Modify facet labels
  facet_wrap(~ Variable, scales = "free_y", labeller = as_labeller(facet_labels)) +
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
p_boxplot_wilcox_outliers_signif



# filtered outliers
p_boxplot_wilcox_rm_outliers <- 
   df_long_narrow_filtered %>% 
  #df_long_narrow %>% 
  mutate(Variable = factor(Variable, 
                           levels = c('prcp', 'tmp',
                                      "distance_edge", "disturbance_severity",
                                      "drought_spei1"
                                     # "mature_dist_severity", 
                                    #  "sum_stems_mature",  
                                     # "residual_mature_trees"
                                    ))) %>% 
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
    size = 0.2#,
    #notch = T
    ) +
  #geom_jitter(size = 0.1, alpha = 0.5, aes(color = adv_delayed)) +
   stat_compare_means(comparisons = comparisons, method = "wilcox.test", 
                     label = "p.signif", 
                     size = 2,
                     label.x = 1.5) +  # Position labels between the groups
  labs(title = "",
       x = "", y = "Vals") +
  theme_classic() +
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
    #panel.grid = element_blank(),          # Optional: Remove grid lines for a cleaner look
    axis.line = element_line(color = "black", linewidth = 0.5)  # Add bottom and left axis lines
  )
#"#A50026" "#DA362A" "#F46D43" "#FDAE61" "#FEE08B" "#D9EF8B" "#A6D96A" "#66BD63" 
# Save the combined plot (optional)
# Save the plot ensuring text sizes are preserved
p_boxplot_wilcox_rm_outliers


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
    min_distance = min(disturbance_severity, na.rm = TRUE),
    max_distance = max(disturbance_severity, na.rm = TRUE),
    mean_distance = mean(disturbance_severity, na.rm = TRUE),
    median_distance = median(disturbance_severity, na.rm = TRUE),
    IQR_distance = IQR(disturbance_severity, na.rm = TRUE),
    sd_distance = sd(disturbance_severity, na.rm = TRUE),
    n = n()
  )


# TEST 

library(ggpubr)
library(dplyr)

# Manually calculate p-values for each comparison
comparison_results <- lapply(comparisons, function(comp) {
  test <- wilcox.test(
    Value ~ adv_delayed, 
    data = df_long_narrow_filtered %>% filter(adv_delayed %in% comp),
    paired = FALSE
  )
  data.frame(
    group1 = comp[1],
    group2 = comp[2],
    p.value = test$p.value
  )
})

# Combine results into a single data frame
comparison_results <- bind_rows(comparison_results)

# Filter for significant results only (p < 0.05)
significant_comparisons <- comparison_results %>%
  dplyr::filter(p.value < 0.05) %>%
  select(group1, group2) %>%
  as.list()

# Create the boxplot
p_boxplot_wilcox_rm_outliers <- 
  df_long_narrow_filtered %>% 
  mutate(Variable = factor(Variable, 
                           levels = c('prcp', 'tmp', "distance_edge", "disturbance_severity",
                                      "mature_dist_severity", 
                                      "sum_stems_mature",  
                                      "residual_mature_trees"))) %>% 
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
    size = 0.2
  ) +
  stat_compare_means(
    #comparisons = significant_comparisons,  # Use only significant comparisons
    method = "wilcox.test", 
    label = "p.format", 
    size = 2,
    label.x = 1.5
  ) +
  labs(title = "",
       x = "", y = "Vals") +
  theme_classic() +
  theme(
    legend.position = 'none',
    text = element_text(size = 8),         # Set all text to size 8
    axis.text = element_text(size = 8),    # Axis tick labels
    axis.title = element_text(size = 8),   # Axis titles
    strip.text = element_text(size = 8),   # Facet labels
    legend.text = element_text(size = 8),  # Legend text
    plot.title = element_text(size = 8),   # Plot title
    strip.background = element_blank(),    # Remove the box around facet names
    strip.placement = "outside",           # Move facet label outside the plot area
    axis.line = element_line(color = "black")  # Add bottom and left axis lines
  )

p_boxplot_wilcox_rm_outliers



# END
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



# Save data ------------------------------------------------
fwrite(df_delayed_advanced, 'outTable/df_delayed_advanced.csv')
fwrite(df_stem_regeneration2, 'outTable/df_stem_regeneration2.csv')


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
library("ggalluvial")


# df_species_summary 
df_species_summary <- 
  df_compare_future_species %>%
  filter(current == 1) %>%  # Keep only species that exist in current conditions
  pivot_longer(cols = c(current, rcp26, rcp45, rcp85), 
               names_to = "scenario", 
               values_to = "presence") %>%
  filter(presence == 1) %>%  # Keep only species that persist
  group_by(acc, scenario) %>%
  summarise(n_sites = n(), .groups = "drop") %>%
  pivot_wider(names_from = scenario, values_from = n_sites, values_fill = 0)  


#df_sankey_test <- 
  df_compare_future_species %>% 
    dplyr::filter(current == 1) %>%
  pivot_longer(c(current, rcp26, rcp45, rcp85), names_to = 'scenario', values_to = 'presence') %>% 
  mutate(scenario = factor(scenario),
         presence = factor(presence)) %>% 
    head()
     
  
  df_freq <- df_compare_future_species %>%
    dplyr::filter(current == 1) %>%
    pivot_longer(c(current, rcp26, rcp45, rcp85), names_to = 'scenario', values_to = 'presence') %>%
    mutate(scenario = factor(scenario),
           presence = factor(presence)) %>%
    group_by(acc, scenario, presence) %>%
    summarise(freq = n(), .groups = "drop")  
  
  
  # Calculate frequency of occurrence grouped by species, scenario, and presence
  df_freq <- df_compare_future_species %>%
    dplyr::filter(current == 1) %>%
    pivot_longer(cols = c(current, rcp26, rcp45, rcp85), names_to = 'scenario', values_to = 'presence') %>%
    mutate(scenario = factor(scenario, levels = c("current", "rcp26", "rcp45", "rcp85")),
           presence = factor(presence)) %>%
    group_by(acc, scenario, presence) %>%
    summarise(freq = n(), .groups = "drop")
  
# try a simple test: still not working, needs further simlification ! 
df_test <- data.frame(site = c(1,1,1,2,2,1,2,2,2,2),
                      site_ch = c("a","a","a","b","b","a","b","b","b","b"),
                      acc = c("piab",
                              "abies",
                              "fasy",
                              "piab",
                              "fasy",
                              "frax",
                              "piab",
                              "fasy",
                              "frax",
                              "abies"
                              ),
                      scenario = rep(c("current", 'rcp25'), each = 5))  %>% 
  mutate(acc = factor(acc),
         site_ch = factor(site_ch),
         scenario = factor(scenario))

# Ensure each species appears across both scenarios per site
df_alluvial <- df_test %>%
  count(site, scenario, acc) %>%  # Count occurrences of each species per site per scenario
  complete(site, scenario, acc, fill = list(n = 0)) %>%  # Ensure all combinations exist
  rename(freq = n) %>%  # Rename count column for clarity
  arrange(site, acc, scenario) %>% 
  dplyr::filter(freq !=0)

# Ensure scenario is a factor in correct order
df_alluvial <- df_alluvial %>%
  mutate(scenario = factor(scenario, levels = c("current", "rcp25")),
         acc = factor(acc),
         site = factor(site))

# Create an alluvial plot
ggplot(df_alluvial,
       aes(x = scenario, stratum = acc, alluvium = site, y = freq, fill = acc)) +
  geom_alluvium(alpha = 0.6, aes(fill = acc)) +  # Connects species across scenarios
  geom_stratum() +  # Adds category blocks
  theme_minimal() +
  labs(title = "Species Presence Across Scenarios",
       x = "Scenario",
       y = "Frequency",
       fill = "Species") +
  scale_fill_brewer(palette = "Set2")





  
  # Ensure scenario is a factor in the correct order
  df_freq <- df_freq %>%
    mutate(scenario = factor(scenario, levels = c("current", "rcp26", "rcp45", "rcp85")),
           presence = factor(presence),
           acc      = factor(acc)) %>% 
    drop_na()
  
  # Create an alluvial plot with scenarios on the x-axis
  ggplot(df_freq,
         aes(x = scenario, stratum = presence, alluvium = acc, y = freq, fill = presence)) +
    #geom_alluvium(aes(fill = presence), alpha = 0.6)# +  # Use geom_alluvium instead of geom_flow
    geom_stratum()# +
    theme_minimal() +
    labs(title = "Species Presence Across Climate Scenarios",
         x = "Scenario",
         y = "Frequency",
         fill = "Presence")
  
head(df_sankey_test)

ggplot(df_sankey_test,
       aes(x = scenario, stratum = presence, alluvium = acc,
           fill = presence, label = presence)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("student curricula across several semesters")







data(majors)
majors$curriculum <- as.factor(majors$curriculum)
ggplot(majors,
       aes(x = semester, stratum = curriculum, alluvium = student,
           fill = curriculum, label = curriculum)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("student curricula across several semesters")




df_sankey <- df_species_summary %>%
  pivot_longer(cols = c(rcp26, rcp45, rcp85), 
               names_to = "scenario", 
               values_to = "n_sites") %>%
  filter(n_sites > 0)  # Remove zero values (no persistence)

# Convert to a format suitable for ggalluvial (Sankey)
df_sankey_alluvial <- df_sankey %>%
  rename(source = current, target = scenario) %>%
  group_by(acc, source, target) %>%
  summarise(flow = sum(n_sites), .groups = "drop")  # Summarize species flow

# Create Sankey plot using ggalluvial
ggplot(df_sankey_alluvial, aes(axis1 = source, axis2 = target, y = flow)) +
  geom_alluvium(aes(fill = acc)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Current", "rcp26", "rcp45", "rcp85")) +
  labs(title = "Species Persistence Sankey Diagram",
       x = "Scenario",
       y = "Number of Sites") +
  theme_minimal()






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

