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

# Read data -----------------------------------------------------------------------

# get vegetation data
load("outData/veg.Rdata")

# final tables on site level
df_fin <- fread('outData/indicators_for_cluster_analysis.csv')

## read coordinates: to add XY coordinates to final model ---------------------
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




## Exploratory analysis ---------------------

### explore general trends of vegetation along management gradient: -------------------------------

# Define the function
create_xy_plot_loess <- function(data, x_var, y_var) {
  data %>%
    ggplot(aes_string(x = x_var, y = y_var)) +  # Use aes_string to allow for dynamic variable names
    geom_jitter(alpha = 0.5) +
    geom_smooth(method = 'loess') +
    theme_bw() +
    theme(aspect.ratio = 1)
}

# Example usage
# Replace 'management_intensity' and 'rIVI' with any desired variables
p1 <- create_xy_plot_loess(df_fin, "management_intensity", "rIVI")
p2 <- create_xy_plot_loess(df_fin, "management_intensity", "richness")
p3 <- create_xy_plot_loess(df_fin, "management_intensity", "stem_density")
p4 <- create_xy_plot_loess(df_fin, "management_intensity", "n_vertical")

# Print the plots (or use in a grid with gridExtra or patchwork)
print(p1)
print(p2)
print(p3)
print(p4)

# add variables:

#library(patchwork)  # Optional: for arranging multiple plots together


# List of variables to plot against 'stem_density'
variables <- c("tmp", "prcp", "tmp_z", "prcp_z", "spei1", "spei3", "spei6", "spei12", "spei24")

# Generate plots and store them in a list
plots <- lapply(variables, function(var) {
  create_xy_plot_loess(df_fin, var, "stem_density")
})

# Optional: Name each plot in the list for easier reference
names(plots) <- paste0("p", 5:(4 + length(variables)))

# Print or save each plot
for (i in seq_along(plots)) {
  print(plots[[i]])  # This will print each plot to the R graphics device
}

# Optional: Combine all plots into a single figure (if needed)
# e.g., using patchwork to arrange all plots together
combined_plot <- ggarrange(plotlist = plots, ncol = 3, nrow = 3)  # Adjust ncol to set the number of columns

# Display the combined plot
(combined_plot)

# Save the combined plot (optional)
ggsave("outFigs/combined_stem_density_plots.png", combined_plot, width = 10, height = 8, dpi = 300)



#View(df_fin)
# get summary for the scatter plot, one point is one country & management
df_summary <- df_fin %>%
  # group_by(manag) %>% # country, 
  summarize(
    rich_med = median(rIVI, na.rm = TRUE),
    rich_sd = sd(rIVI, na.rm = TRUE),
    rich_25 = quantile(rIVI, 0.25, na.rm = TRUE), # rich_mean -rich_sd, #
    rich_75 = quantile(rIVI, 0.75, na.rm = TRUE), # rich_mean +rich_sd,
    dens_med = median(stem_density   , na.rm = TRUE),
    dens_sd   = sd(stem_density  , na.rm = TRUE),
    dens_25 = quantile(stem_density  , 0.25, na.rm = TRUE), # dens_mean - dens_sd, #
    dens_75 = quantile(stem_density  , 0.75, na.rm = TRUE), # dens_mean + dens_sd, #
    .groups = 'drop'
  )

(df_summary)








# 2. species composition per stems: --------------

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



# ad indicators for clim cluster
stem_dens_species_long_cluster <- stem_dens_species_long_cluster %>% 
  left_join(clim_cluster_indicator, by = c('cluster' = 'site')) #%>% 
  

# get overall number for the species compositions:
stem_dens_species_long_cluster %>% 
  group_by(Species) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup(.) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share       = sum_stems/sum_vegType*100)#%>% 
  #group_by(clim_class) %>% 
  
 # arrange(desc(share)) #%>% 
#View()


# get density plot per vertical class and species

# Summarize the data by VegType and Species
df_sum_stems_vertical <- stem_dens_species_long_cluster %>%
  group_by(cluster,VegType, Species) %>% #cluster, 
  summarize(sum_stems = sum(stem_density, na.rm = TRUE)) %>%
  ungroup() %>% 
  dplyr::filter(sum_stems > 0) #%>%  
 # arrange(sum_stems) 


# Function to cap values at 1st and 99th percentiles
cap_percentiles <- function(x, lower = 0.01, upper = 0.99) {
  q <- quantile(x, probs = c(lower, upper), na.rm = TRUE)
  x <- ifelse(x < q[1], q[1], x)
  x <- ifelse(x > q[2], q[2], x)
  return(x)
}

# Cap the sum_stems data at 1st and 99th percentiles for each Species and VegType
df_sum_stems_vertical_capped <- df_sum_stems_vertical %>%
  group_by(Species, VegType) %>%
  mutate(sum_stems_capped = cap_percentiles(sum_stems)) %>%
  ungroup()

  # plot just meadian and IQR
  # test plotting START -----------------
  # Create the density plot using ggplot2

# find dominant/prevailing species: (this is now an example!)

dominant_species5 <- c('piab', 'fasy', 'pisy', 'besp')
  df_test <- stem_dens_species_long_cluster %>% 
    dplyr::filter(stem_density > 0) %>%  
    dplyr::filter(Species %in% dominant_species5) #%>% 
  
df_test %>% 
    ggplot(aes(x = VegType, y = stem_density, fill = VegType)) +
    geom_boxplot(outlier.shape = NA) +  # Hide outliers
  geom_jitter(size = 0.5, alpha = 0.2) +
    coord_flip() +  # Flip the coordinates
    ylim(0, 5000) +  # Zoom in on the stem_density range
    facet_grid(Species ~ .) +  # Facet by Species
   theme_classic() +
  theme(legend.position = 'none')
  
df_test %>% 
  ggplot(aes(x = VegType, y = stem_density, fill = VegType)) +
  geom_violin(trim = T) +  # trim extreme values
  geom_jitter(size = 0.5, alpha = 0.2) +
  coord_flip() +  # Flip the coordinates
  ylim(0, 5000) +  # Zoom in on the stem_density range
  facet_grid(Species ~ .) +  # Facet by Species
  theme_classic() +
  theme(legend.position = 'none')

  
  
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
  
  
# Species composition: vertical class  ------------------------------------
# Summarize the data by VegType and Species
mean_stems_vertical <- stem_dens_species_long %>%
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
# Summarize the total stem density per species for each climate class
species_composition <- stem_dens_species_long %>%
  group_by(clim_class, Species) %>%
  summarize(sum_stems = sum(stem_density, na.rm = TRUE)) %>% 
  ungroup() 

# Calculate the total stem density per climate class and the share of each species
species_composition <- species_composition %>%
  group_by(clim_class) %>%
  mutate(total_stems_clim_class = sum(sum_stems),  # Total stem density in each climate class
         share = (sum_stems / total_stems_clim_class) * 100) %>%  # Calculate percentage share
  ungroup()


#sum of species per fiel wrork: from 37 tree species
# Calculate species richness per clim_class
species_richness <- top_species_per_clim_class %>%
  group_by(clim_class) %>%
  summarise(species_richness = n_distinct(Species))

stem_dens_species_long %>% 
  ungroup() %>% 
  dplyr::filter(stem_density > 0) %>% 
  summarise(species_richness = n_distinct(Species))
 


# Use different gradient depepnding fof teh seral stage: 

# Find the top 5 species per climate class based on share
top_species_per_clim_class <- species_composition %>%
  group_by(clim_class) %>%
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
species_composition <- stem_dens_species_long %>%
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
top_species_per_clim_class <- species_composition %>%
  group_by( country) %>%
  arrange(desc(share)) %>%  # Sort species by their share within each climate class
  #slice_head(n = 6) #%>%  # Select the top 5 species per country
  dplyr::filter(share > 5)# %>%  # select species with share > 5%
  #left_join(df_seral_species, by = join_by(Species))  #%>% 
#arrange(seral_type) %>% 
#mutate(seral_type=factor(seral_type) )

# Ensure species are arranged by seral type and within each climate class
top_species_per_clim_class <- top_species_per_clim_class %>%
  arrange(country, Species)  # First arrange by seral type and then by Species alphabetically

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

stem_dens_species_long %>% 
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

## Get presence absence data for indiviual layers --------------------------

# Group by cluster and VegType, check if stem_density > 0
vert_class_presence_absence <- stem_dens_species_long %>%
  group_by(cluster, VegType) %>%
  summarise(has_stem_density = sum(stem_density > 0, na.rm = TRUE) > 0) %>%
  ungroup() %>%
  # Filter only where there is stem density present
  dplyr::filter(has_stem_density) %>%
  dplyr::select(cluster, VegType)

# Now find clusters where no VegType has stem_density > 0
# Find clusters where no vegType has stem_density
all_clusters <- stem_dens_species_long %>%
  ungroup() %>% 
  dplyr::select(cluster, VegType) %>% 
  distinct(cluster)

clusters_with_stem_density <- presence_absence %>%
  ungroup() %>% 
  dplyr::select(cluster, VegType) %>% 
  distinct(cluster)

clusters_without_stem_density <- all_clusters %>%
  anti_join(clusters_with_stem_density, by = "cluster") %>%
  mutate(VegType = "0")

# Combine clusters with and without stem density
final_result <- bind_rows(vert_class_presence_absence, clusters_without_stem_density) %>%
  arrange(cluster, VegType)

# Display final result
final_result

# reclassify naming convention

final_result <- final_result %>%
  mutate(VegType = case_when(
    VegType == "Mature" ~ "m",
    VegType == "Saplings" ~ "s",
    VegType == "Juveniles" ~ "j",
    TRUE ~ VegType # In case there are other values
  )) %>% 
  left_join(clim_cluster_indicator, by = c('cluster' = 'site'))  

# Group by cluster, concatenate VegType into a single string for each cluster
vert_vegtype_classes <- final_result %>%
  group_by(cluster, clim_class) %>%
  summarise(veg_class = paste(sort(unique(VegType)), collapse = "")) %>%
  ungroup()

# Count the occurrences of each veg_class within each clim_class
vert_vegtype_class_count <- vert_vegtype_classes %>%
  group_by(clim_class, veg_class) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  arrange(clim_class, veg_class)

# Calculate the total counts per clim_class
vert_vegtype_class_percent <- vert_vegtype_class_count %>%
  group_by(clim_class) %>%
  mutate(total_count = sum(count)) %>% 
  ungroup() %>%
  # Calculate the percentage for each veg_class within each clim_class
  mutate(percent = round((count / total_count) * 100, 1))


# order the classes:

# Set up the correct factor levels for veg_class based on your preferred order
vert_vegtype_class_percent$veg_class <- factor(vert_vegtype_class_percent$veg_class, 
                                               levels = c("0", "s", "js", "j", "jm",  "jms", "ms", "m"))

# Define a green color palette using colorRampPalette for 's' to 'jms' classes
green_palette <- colorRampPalette(c("lightgreen", "darkgreen"))(7)  # 7 shades of green

# Combine the red color for '0' with the generated green palette
color_palette <- c("0" = "red", setNames(green_palette, c("s", "js", "j", "jm",  "jms", "ms", "m")))



# Create a stacked bar plot
p_vertical_classes <- ggplot(vert_vegtype_class_percent, aes(x = clim_class, y = percent, fill = veg_class)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar plot
  geom_text(aes(label = #paste0(count, " (", percent, "%)")
                  ifelse(percent>=2, percent, '')
  ),
  position = position_stack(vjust = 0.5),  # Labels inside the bars
  size = 3, color = "black") +  # Adjust text size and color
  labs(x = "", y = "Percentage", 
       fill = "Vertical\nclass",
       title = "") +
  scale_fill_manual(values = color_palette) +  # Apply the color palette with green gradient and red
  theme_classic2() +  # Use a clean theme
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    plot.title = element_text(hjust = 0.5)  # Center the title
  )



ggsave(filename = 'outFigs/fig_vert_classes_share_clim_clust.png', 
       plot = p_vertical_classes, 
       width = 5, height = 4, dpi = 300, bg = 'white')



# 


#### test vertical layers - OLD -----------------------------------------------------

# Get again individual layers:
df_vert_full <- 
  stem_dens_species_long %>% 
  left_join(clim_cluster_indicator, by = c('cluster' = 'site')) %>% 
  dplyr::filter(stem_density>0) %>% 
  ungroup(.) %>% 
  dplyr::select(cluster, VegType) %>%
  group_by(cluster) %>%
  distinct(.) %>% 
  mutate(n_layers = n()) #%>% 



# Calculate frequencies of individual layers and combinations
layer_frequencies <- 
  df_vert_full %>%
  group_by(cluster) %>% 
  right_join(df_stems) %>% # to account for teh empty ones
  mutate(VegType = case_when(
    is.na(VegType) ~ 'NO',   # Replace NA with 'NO'
    TRUE ~ VegType           # Keep existing values for non-NA cases
  )) %>% 
  summarise(layers = paste(sort(unique(VegType)), collapse = "-")) %>%
  ungroup() %>%
  count(layers)

layer_frequencies$prop <- layer_frequencies$n/n_clusters*100
(layer_frequencies)

# plot how often each category occurs
windows()
layer_frequencies %>% 
  ggplot(aes(x = reorder(layers, -n), y = n)) +
  geom_bar(stat = "identity") #, position = "fill") #+
#geom_text(aes(label = paste0(round(scales::percent, 1), "%"), 
#              y = cumsum(n) - n/2), 
#          position = position_fill(), 
#          color = "white", 
#          size = 3) +
ylab("Percentage (%)") +
  theme_classic() +
  labs(fill = "Layers")


layer_frequencies %>% 
  ggplot(aes(fill = layers     , 
             area = prop      ,
             label = paste(layers, "\n", round(prop,2), "%"))) +
  geom_treemap() +
  geom_treemap_text(colour ="white", place = "centre")# +




layer_reg_frequencies <- 
  df_vert_full %>%
  dplyr::filter(VegType != 'Mature') %>% 
  group_by(cluster) %>% 
  right_join(df_stems) %>% # to account for teh empty ones
  mutate(VegType = case_when(
    is.na(VegType) ~ 'NO',   # Replace NA with 'NO'
    TRUE ~ VegType           # Keep existing values for non-NA cases
  )) %>% 
  summarise(layers = paste(sort(unique(VegType)), collapse = "-")) %>%
  ungroup() %>%
  count(layers)


# View the result
print(layer_reg_frequencies)


###### chi sqare goodnes of fit: -------------------
# expecting that categories are equally distributed: divide total observations by number of groups
# Total observations
total_observations <- sum(layer_reg_frequencies$n)

# Expected frequencies if all layers were equally likely
expected_frequencies <- rep(total_observations / nrow(layer_reg_frequencies), 
                            nrow(layer_reg_frequencies))

# Perform the chi-square test for goodness of fit
chi_square_test_result <- chisq.test(x = layer_reg_frequencies$n, 
                                     p = expected_frequencies / total_observations)

# Print the result
print(chi_square_test_result)


# compare residuals to understand the differences between groups:
# Assuming chi_square_test_result is the result of your chi-square test
# and layer_reg_frequencies contains the observed frequencies for your categories

# Calculate expected frequencies
expected_frequencies <- chi_square_test_result$expected

# Calculate standardized residuals
standardized_residuals <- (layer_reg_frequencies$n - expected_frequencies) / sqrt(expected_frequencies)

# Associate residuals with group names
residuals_table <- data.frame(Group = layer_reg_frequencies$layers, 
                              StandardizedResiduals = standardized_residuals)

# Print the residuals table, sorted by the absolute value of residuals
residuals_table <- residuals_table[order(abs(residuals_table$StandardizedResiduals), decreasing = TRUE),]
print(residuals_table)











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


# visualize vertical classes by UpSEt plot

# Upset plot -----------------------------------

library(UpSetR)

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



# 



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
# split table in two: drivers for advanced (> 1000 stems/ha)
#                     drivers for delayed regeneration (<50 stems/ha)  

# Prep fnal table --------------

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

# make two regeneration classes: advnaced vs advanced
df_fin <- df_fin %>% 
  mutate(adv_delayed = ifelse(stem_regeneration <= 50, "Delayed", 
                              ifelse(sum_stems_juvenile >= 1000, "Advanced", NA))) %>% 
  mutate(delayed = ifelse(stem_regeneration <= 50, 1, 0),
         advanced = ifelse(sum_stems_juvenile >=  1000, 1, 0))





## Variables selection: Corelations  ------------------------------------------------------

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
library(Boruta)
library(randomForest)
 
df_fin$log_distance_edge        <- log(df_fin$distance_edge + 1)
df_fin$log_disturbance_severity <- log(df_fin$disturbance_severity + 0.001)
df_fin$log_depth_extract        <- log(df_fin$depth_extract + 1) # in cm
df_fin$salvage_intensity        <- df_fin$salvage_intensity + 0.01 # increase a bit


# occurence or risk of teh binary outcome
# Create the binary outcome variable for regeneration status
df_fin$delayed <- ifelse(df_fin$stem_regeneration < 50, 1,0)
df_fin$advanced <- ifelse(df_fin$stem_regeneration > 1000, 1,0)
df_fin$advanced_juv <- ifelse(df_fin$sum_stems_juvenile > 1000, 1,0)





hist(df_fin$log_depth_extract)
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
                    "drought_spei12", 
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



# Subset the dataset to include only predictor and dependent variables
df_for_boruta <- df_fin[, c(dependent_vars, predictor_vars_sub), with=FALSE]

# check for correlations and remote the highly correlated predictors
cor_matrix <- cor(df_for_boruta[, predictor_vars_sub, with=FALSE], use="pairwise.complete.obs")
# View correlation matrix to find high correlations
round(cor_matrix, 2)

# remove high correlated predictor:
# Set the upper triangle and the diagonal of the correlation matrix to NA
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA

# Find pairs of variables with correlations above the threshold (e.g., 0.75)
high_cor_pairs <- which(abs(cor_matrix) > 0.7, arr.ind = TRUE)

# Display highly correlated pairs
print(high_cor_pairs)

# Identify the names of predictors to be removed (keep the first occurrence)
predictors_to_remove <- unique(rownames(high_cor_pairs))

predictor_vars_cleaned <- predictor_vars_sub[!(predictor_vars_sub %in% predictors_to_remove)]

# Display cleaned predictor variables
print(colnames(predictor_vars_cleaned))



# TEST LOOP BORUTA ---------------------------------------

# Define the function that performs Boruta analysis and returns the merged result
run_boruta_analysis <- function(dep_var, df_for_boruta, predictor_vars_sub) {
  
  # Run Boruta analysis
  boruta_results <- Boruta(as.formula(paste(dep_var, "~ .")),
                           data = df_for_boruta[, c(dep_var, predictor_vars_sub), with = FALSE], 
                           doTrace = 2)
  
  # Extract final decisions for each variable
  decisions_df <- as.data.frame(boruta_results$finalDecision)
  decisions_df$Dependent <- dep_var  # Create new column for the dependent variable
  
  # Add row names as a new column called "Variable"
  decisions_df$Variable <- rownames(decisions_df)
  
  # Rename the decision column to "Decision"
  colnames(decisions_df) <- c("Decision", "Dependent", "Variable")
  
  # Extract importance history and remove shadow variables
  importance_history <- boruta_results$ImpHistory[, !grepl("shadow", colnames(boruta_results$ImpHistory))]
  
  # Convert the importance history to a data frame
  importance_df <- as.data.frame(importance_history)
  
  # Pivot the data to long format
  importance_long <- importance_df %>%
    pivot_longer(cols = everything(),  # Select all columns
                 names_to = "Variable", 
                 values_to = "Importance")
  
  # Merge the importance scores with the Boruta decisions
  result_df <- importance_long %>%
    left_join(decisions_df, by = "Variable")
  
  return(result_df)
}

# Apply the function to each dependent variable using lapply or a loop
out_list <- lapply(dependent_vars, function(dep_var) {
  print(dep_var)
  run_boruta_analysis(dep_var, df_fin, predictor_vars_cleaned)
})

# Combine the results into a single dataframe
out_df <- bind_rows(out_list)

# View the final output
head(out_df)

# Order variables by median importance, from lowest to highest
out_df <- out_df %>%
  group_by(Variable) %>%
  mutate(median_importance = median(Importance, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(median_importance)

# Create a factor for ordered variables by median importance
out_df$Variable <- factor(out_df$Variable, levels = unique(out_df$Variable))

p_select_vars_boruta <- out_df %>% 
  dplyr::filter(!grepl('sum_', Dependent)) %>% 
  ggplot(aes(x = Variable, y = Importance, fill = Decision)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Confirmed" = "green", "Rejected" = "red", "Tentative" = "orange")) +
  labs(title = "Boxplot of Importance by Variables",
       x = "Variables", 
       y = "Importance") +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') + # Rotate x-axis labels for better readability
facet_wrap(Dependent ~.)

ggsave('outFigs/fig_select_vars_boruta.png', plot = p_select_vars_boruta, width =7,height = 7)
# get sumary table:

out_df_summary <- out_df %>% 
       dplyr::filter(Importance != -Inf) %>% 
  group_by(Dependent, Variable, Decision ) %>% 
  summarise(median = median(Importance , na.rm = T))# %>% 
  #group_by(Dependent) %>% 
  #slice_max(order_by = median, n = 5) %>%
  #print(n = 40)
# 

View(out_df_summary)




## Models: prepare fin tables ---------------------------------------------------------------------------------

#
# make median +IQR plots of variables for quick view


# Subsetting the data for different models
# Variables to be used for each dependent variable
predictors_delayed <- c("rIVI", 
                        # "n_vertical",  # remove as it has only 1/0 outcome
                        # "richness",    # remove as it has only 1/0 outcome
                        "clay_extract",
                        "sum_stems_mature", 
                        "prcp",
                        "spei12",
                        "disturbance_severity",
                        "distance_edge",
                        #"drought_spei12",
                        "tmp",
                        "av.nitro",
                        "depth_extract",
                        "management_intensity",
                        "country_pooled","clim_grid", "clim_class", "x", "y")


predictors_advanced <- c("richness", 
                         "rIVI", 
                         "n_vertical", 
                         "sum_stems_mature", 
                         "spei12",
                        # "drought_spei12", 
                         "country_pooled",
                         "prcp",
                         "tmp",
                         "clay_extract",
                         "av.nitro",
                         "depth_extract",
                         "management_intensity",
                         "disturbance_severity",
                        "distance_edge",
                         "clim_grid",  
                         "clim_class", 
                         "x", "y")
predictors_stem_regeneration <- c("richness",
                                  "n_vertical", 
                                  "rIVI", 
                                  "prcp", 
                                  #"drought_spei12",
                                  "spei12",
                                  "tmp",
                                  "av.nitro",
                                  "clay_extract",
                                  "depth_extract",
                                  "management_intensity",
                                  "disturbance_severity",
                                  "distance_edge",
                                  "country_pooled", "clim_grid", "clim_class", "x", "y")

# Subset the data
df_delayed <- df_fin %>% dplyr::select(all_of(c("delayed", predictors_delayed)))
df_advanced <- df_fin %>% dplyr::select(all_of(c("advanced", predictors_advanced)))
df_stem_regeneration <- df_fin %>% dplyr::select(all_of(c("stem_regeneration", predictors_stem_regeneration)))



# test drivers: simplify the analysis:
# Subset the data
df_delayed2 <- df_fin %>% 
  dplyr::select(all_of(c("delayed", predictor_vars_sub, "management_intensity",
                                                "country_pooled", "clim_grid", "clim_class", "x", "y", "tmp_c", "prcp_c")))
df_advanced2 <- df_fin %>% 
  dplyr::select(all_of(c("advanced", predictor_vars_sub, "management_intensity",
                                                 "country_pooled", "clim_grid", "clim_class", "x", "y", "tmp_c", "prcp_c")))

#df_advanced_juv2 <- df_fin %>% 
#  dplyr::select(all_of(c("advanced_juv", predictor_vars_sub, "management_intensity",
 #                        "country_pooled", "clim_grid", "clim_class", "x", "y")))

df_stem_regeneration2 <- df_fin %>% 
  dplyr::select(all_of(c("stem_regeneration", predictor_vars_sub,"management_intensity",
                                                          "country_pooled", "clim_grid", "clim_class", "x", "y")))


table(df_fin$delayed)
table(df_fin$advanced)




# Drivers: ---------------------------------------------------------------------


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





##  Drivers delayed  ---------------------------

# Check the structure of the data
str(df_delayed2)

# Check for missing values
sum(is.na(df_delayed2))

# Inspect the distribution of each predictor
hist(df_delayed2$prcp, main="Histogram of Precipitation (prcp)", xlab="prcp")
hist(df_delayed2$tmp, main="Histogram of Temperature (tmp)", xlab="tmp")

# Check collinearity using the VIF
vif_check <- performance::check_collinearity(df_delayed2)
print(vif_check)

# # END
# # Adding a small constant to handle zero-inflation
# df_delayed2$richness_adj <- df_delayed2$richness + 0.01
# df_delayed2$n_vertical_adj <- df_delayed2$n_vertical + 0.01
# df_delayed2$sum_stems_mature_adj <- df_delayed2$sum_stems_mature + 0.01


library(corrplot)
numeric_vars <- df_delayed2[, sapply(df_delayed2, is.numeric)]
cor_matrix <- cor(na.omit(numeric_vars))
corrplot(cor_matrix)

ggplot(df_delayed2, aes(x = richness)) + geom_histogram(bins = 20) + facet_wrap(~ delayed)

delayed_basic_gam <- gam(delayed ~ s(prcp, k = 10) + s(tmp, k = 10) + s(spei12, k = 10) + 
                   s(distance_edge) +
                   s(depth_extract) +
                   s(disturbance_severity) +
                   s(clay_extract) +
                   s(av.nitro) +
                   #s(richness_adj, k = 3) + 
                   #s(n_vertical_adj, k = 3) + 
                   #s(sum_stems_mature_adj, k = 5) + 
                   s(management_intensity, k = 5) +
                   s(country_pooled, bs = 're') # + clim_grid + clim_class
                 , 
                 family = binomial(link = "logit"), data = df_delayed2, method = "REML")
summary(basic_gam)
gam.check(basic_gam)
appraise(basic_gam)
windows()
plot.gam(basic_gam, page = 1)


delayed_interaction_model_1 <- gam(delayed ~ 
                             s(prcp, k = 10) + 
                             s(tmp, k = 10) + 
                             s(spei12, k = 10) +
                             s(distance_edge) + 
                             s(depth_extract) +
                             s(disturbance_severity) +
                             s(clay_extract) + 
                             s(av.nitro) +
                             s(management_intensity, k = 5) +
                             s(country_pooled, bs = "re") +
                             ti(prcp, tmp, k = 5), 
                           family = binomial(link = "logit"), 
                           data = df_delayed2, 
                           method = "REML")

summary(interaction_model_1)

delayed_interaction_model_2 <- gam(delayed ~ 
                             s(prcp, k = 10) + 
                             s(tmp, k = 10) + 
                             s(spei12, k = 10) +
                             s(distance_edge) + 
                             s(depth_extract) +
                             s(disturbance_severity) +
                             s(clay_extract) + 
                             s(av.nitro, k = 15) +
                             s(management_intensity, k = 5) +
                             s(country_pooled, bs = "re") +
                             s(clim_grid, bs = "re") +
                             ti(management_intensity, disturbance_severity) +
                             #ti(tmp, spei12, k = 5) +
                             ti(prcp, tmp, k = 5), 
                           family = binomial(link = "logit"), 
                           data = df_delayed2, 
                           method = "REML")

delayed_interaction_model_2_upd <- gam(delayed ~ 
                                     s(prcp, k = 15) + 
                                     s(tmp, k = 10) + 
                                     s(spei12, k = 10) +
                                     s(distance_edge) + 
                                     s(depth_extract) +
                                     s(disturbance_severity) +
                                     s(clay_extract) + 
                                     s(av.nitro, k = 15) +
                                       s(management_intensity, by = country_pooled, bs = "re") +
                                     #s(management_intensity, k = 5) +
                                     s(country_pooled, bs = "re") +
                                     s(clim_grid, bs = "re") +
                                     ti(management_intensity, disturbance_severity) +
                                     #ti(tmp, spei12, k = 5) +
                                     ti(prcp, tmp, k = 5) +
                                       s(x, y)  , 
                                   family = binomial(link = "logit"), 
                                   data = df_delayed2, 
                                   method = "REML")

# remove interaction between management and disturbance intensity - difficult to explain
delayed_interaction_model_2_upd2 <- gam(delayed ~ 
                                         s(prcp, k = 15) + 
                                         s(tmp, k = 10) + 
                                         s(spei12, k = 10) +
                                         s(distance_edge) + 
                                         s(depth_extract) +
                                         s(disturbance_severity) +
                                         s(clay_extract) + 
                                         s(av.nitro, k = 15) +
                                         s(management_intensity, by = country_pooled, bs = "re") +
                                         #s(management_intensity, k = 5) +
                                         s(country_pooled, bs = "re") +
                                         s(clim_grid, bs = "re") +
                                         #ti(management_intensity, disturbance_severity) +
                                         #ti(tmp, spei12, k = 5) +
                                         ti(prcp, tmp, k = 5) +
                                         s(x, y)  , 
                                       family = binomial(link = "logit"), 
                                       data = df_delayed2, 
                                       method = "REML")


# disturbance severity and spei seems intersting: add one more term!
# remove interaction between management and disturbance intensity - difficult to explain
delayed_interaction_model_2_upd3 <- gam(delayed ~ 
                                          s(prcp, k = 15) + 
                                          s(tmp, k = 10) + 
                                          s(spei12, k = 10) +
                                          s(distance_edge) + 
                                          s(depth_extract) +
                                          s(disturbance_severity) +
                                          s(clay_extract) + 
                                          s(av.nitro, k = 15) +
                                          s(management_intensity, by = country_pooled, bs = "re") +
                                          #s(management_intensity, k = 5) +
                                          s(country_pooled, bs = "re") +
                                          s(clim_grid, bs = "re") +
                                          #ti(management_intensity, disturbance_severity) +
                                          ti(disturbance_severity, spei12, k = 5) +
                                          ti(prcp, tmp, k = 5) +
                                          s(x, y)  , 
                                        family = binomial(link = "logit"), 
                                        data = df_delayed2, 
                                        method = "REML")




# no xy
delayed_interaction_model_2_upd_no_xy <- gam(delayed ~ 
                                         s(prcp, k = 10) + 
                                         s(tmp, k = 10) + 
                                         s(spei12, k = 10) +
                                         s(distance_edge) + 
                                         s(depth_extract) +
                                         s(disturbance_severity) +
                                         s(clay_extract) + 
                                         s(av.nitro, k = 15) +
                                         s(management_intensity, by = country_pooled, bs = "re") +
                                         #s(management_intensity, k = 5) +
                                         s(country_pooled, bs = "re") +
                                         s(clim_grid, bs = "re") +
                                         ti(management_intensity, disturbance_severity) +
                                         #ti(tmp, spei12, k = 5) +
                                         ti(prcp, tmp, k = 5) #+
                                         #s(x, y) 
                                         , 
                                       family = binomial(link = "logit"), 
                                       data = df_delayed2, 
                                       method = "REML")

# no xy
delayed_interaction_model_2_upd_no_xy_interct2 <- gam(delayed ~ 
                                               s(prcp, k = 10) + 
                                               s(tmp, k = 10) + 
                                               s(spei12, k = 10) +
                                               s(distance_edge) + 
                                               s(depth_extract) +
                                               s(disturbance_severity) +
                                               s(clay_extract) + 
                                               s(av.nitro, k = 15) +
                                               s(management_intensity, by = country_pooled, bs = "re") +
                                               #s(management_intensity, k = 5) +
                                               s(country_pooled, bs = "re") +
                                               s(clim_grid, bs = "re") +
                                               ti(management_intensity, disturbance_severity) +
                                               ti(tmp, spei12, k = 5) +
                                               ti(prcp, tmp, k = 5) #+
                                             #s(x, y) 
                                             , 
                                             family = binomial(link = "logit"), 
                                             data = df_delayed2, 
                                             method = "REML")



AIC(delayed_interaction_model_2_upd, delayed_interaction_model_2,delayed_interaction_model_2_upd2,
    delayed_interaction_model_2_upd3,
    delayed_interaction_model_2_upd_no_xy,delayed_interaction_model_2_upd_no_xy_interct2)

summary(delayed_interaction_model_2_upd3)
summary(delayed_interaction_model_2)
appraise(delayed_interaction_model_2)
appraise(delayed_interaction_model_2_upd_no_xy)

delayed_interaction_model_3 <- gam(delayed ~ 
                             s(prcp, k = 10) + 
                             s(tmp, k = 10) + 
                             s(spei12, k = 10) +
                             s(distance_edge) + 
                             s(depth_extract) +
                             s(disturbance_severity) +
                             s(clay_extract) + 
                             s(av.nitro) +
                             s(management_intensity, k = 5) +
                             s(country_pooled, bs = "re") +
                             ti(tmp, management_intensity, k = 5), 
                           family = binomial(link = "logit"), 
                           data = df_delayed2, 
                           method = "REML")

summary(delayed_interaction_model_3)
AIC(delayed_basic_gam, delayed_interaction_model_1, delayed_interaction_model_2, delayed_interaction_model_3)
BIC(delayed_basic_gam, delayed_interaction_model_1, delayed_interaction_model_2, delayed_interaction_model_3)

predicted_interaction <- ggpredict(delayed_interaction_model_2, terms = c("prcp", "spei12"))
plot(predicted_interaction)

windows()
plot.gam(fin.m.delayed, page = 1)
appraise(delayed_interaction_model_2)
summary(delayed_interaction_model_2)
gam.check(delayed_interaction_model_2)
k.check(delayed_interaction_model_2)

# see influential points:
influence_data <- influence.gam(interaction_model_2)
summary(influence_data)

# identify observations with high influence
high_influence <- which(influence_data > 0.2)  # Adjust the threshold as needed
df_delayed2[high_influence, ]


draw(interaction_model_2)
summary(interaction_model_2)

# get morans stats to test spatial autocorrelation

summary(delayed_interaction_model_2_upd)
gam.check(delayed_interaction_model_2_upd)
k.check(delayed_interaction_model_2_upd)
appraise(delayed_interaction_model_2_upd)
appraise(delayed_interaction_model_2)


# save up teh best model
fin.m.delayed <- delayed_interaction_model_2_upd

# Test for spatial autocorrelation
model_residuals <- residuals(fin.m.delayed, type = "pearson")

# Load necessary libraries
library(spdep)

# Create coordinates matrix
coords <- cbind(df_delayed2$x, df_delayed2$y)

# Create a spatial neighbors object (e.g., using k-nearest neighbors)
# Adjust k based on the density and distribution of your data points
nb <- knn2nb(knearneigh(coords, k = 5))

# Convert neighbors list to a weights list
listw <- nb2listw(nb, style = "W")
# Perform Moran's I test
moran_test <- moran.test(model_residuals, listw)
print(moran_test)



### Plot: drivers delayed ---------------------------------------------------
summary(fin.m.delayed)

y_lab = 'Risk delayed'

# Predict management intensity effects by country
predicted_mgmt_intensity <- ggpredict(fin.m.delayed, terms = c("management_intensity", "country_pooled"))
# Generate predictions for stem regeneration across different countries
predicted_countries <- ggpredict(fin.m.delayed, terms = "country_pooled")

# Create plots for individual variables
plot_prcp   <- create_plot(fin.m.delayed, term = "prcp", df_delayed2 , "Prcp *", scatter_y = "delayed",line_color = "blue", fill_color = "blue")
plot_tmp    <- create_plot(fin.m.delayed, "tmp", df_delayed2, "Tmp **", scatter_y = "delayed",line_color = "red", fill_color = "red")
#plot_spei12 <- create_plot(fin.m.delayed, "spei12", df_delayed2, "spei12", scatter_y = "delayed",line_color = "green", fill_color = "green")
#plot_disturbance_severity <- create_plot(fin.m.delayed, "disturbance_severity", df_delayed2, "Dist Sev", scatter_y = "delayed", line_color = "purple", fill_color = "purple")
#plot_distance_edge <- create_plot(fin.m.delayed, "distance_edge", df_delayed2, "Distance edge", x_limit = c(0, 600), scatter_y = "delayed", line_color = "grey", fill_color = "grey")
#plot_management_intensity <- create_plot(fin.m.delayed, "management_intensity", df_delayed2, "Mng intend",scatter_y = "delayed", line_color = "grey", fill_color = "grey")


# Create interaction plots
plot_interaction1 <- create_interaction_plot(fin.m.delayed, 
                                             c("prcp", "tmp[8,9,10]"), "prcp & tmp *", df_delayed2) +
  guides(color = guide_legend(title = "tmp"), fill = guide_legend(title = "tmp")) + # Rename legend to 'tmp'
  theme(legend.position.inside = c(0.7, 0.7),
    #legend.position.inside = "top right", # Place the legend in the upper right corner inside the plot
        legend.background = element_rect(fill = "white", color = NA), # Add a white background
        legend.title = element_text(size = 10), # Adjust the title size
        legend.text = element_text(size = 8))   # Adjust the text size

plot_interaction2 <- create_interaction_plot(fin.m.delayed, c("management_intensity", "disturbance_severity[0.5,0.9]"), "manag & dist severity *", 
                                             df_delayed2) +
  guides(color = guide_legend(title = "tmp"), fill = guide_legend(title = "tmp")) + # Rename legend to 'tmp'
  theme(legend.position.inside =  c(0.7, 0.7), # Place the legend in the upper right corner inside the plot
        legend.background = element_rect(fill = "white", color = NA), # Add a white background
        legend.title = element_text(size = 10), # Adjust the title size
        legend.text = element_text(size = 8))   # Adjust the text size

# Plot predicted effects by country
plot_countries <- ggplot(predicted_countries, aes(x = x, y = predicted, fill = x)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = "", x = "Country", y = y_lab) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the effect of management intensity on stem regeneration by country
plot_mgmt_intensity <- ggplot(predicted_mgmt_intensity, aes(x = x, y = predicted, group = group)) +
  geom_line(aes(color = group), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  labs(title = "", 
       x = "Management Intensity", y = y_lab, color = "Country", fill = "Country") +
  theme_minimal() +
  facet_wrap(~group) +
  theme(legend.position = "none")


# !!!
library(cowplot)
#plot_grid(p1, p2, labels = "auto", label_size = 12)

# Combine the individual plots into one figure
combined_plot <- plot_grid(plot_interaction1, plot_interaction2, #combined_interactions,
                           #combined_countries,
                           plot_countries,NULL,
                            plot_prcp, plot_tmp,
                          #plot_spei12, 
                          #plot_disturbance_severity, 
                          #plot_distance_edge,             
                          #plot_management_intensity,
                           ncol = 2, nrow = 3,
                          rel_heights = c(1, 1,1))

combined_plot

ggsave('outFigs/fig_delayed_drivers.png', plot = combined_plot, width =5,height = 5, bg = 'white')
# Compute explained deviance
explained_deviance <- round(100 * summary(fin.m.delayed)$dev.expl, 2)  # Calculate explained deviance in percentage

sjPlot::tab_model(fin.m.delayed, 
                  file = "outTable/full_drivers_delayed.doc",
                  dv.labels = paste0("Explained Deviance: ", explained_deviance, "%"),
                  show.re.var = TRUE, # show random effects
                  show.dev = F)    # show deviance explained


# test
# Fit your model
model <- delayed_interaction_model_2_upd2

# Identify random effects using the model's "smooth" component
smooth_terms <- summary(model)$s.table

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
sjPlot::tab_model(model,
          show.re.var = TRUE,        # Show the variance components
          #show.icc = TRUE,           # Show Intraclass Correlation Coefficient
          #show.dev = TRUE,           # Show deviance
          pred.labels = c("Intercept", pred_labels), # Replace smooth term labels
          dv.labels = paste0("Explained Deviance: ", round(100 * summary(model)$dev.expl, 2), "%"), 
          file = "outTable/full_drivers_delayed.doc")


# end



## drivers advanced regen ---------------------------------------------------------

# Check the structure of the data
str(df_advanced2)

# Check for missing values
sum(is.na(df_advanced2))

# Inspect the distribution of each predictor
hist(df_advanced2$prcp, main="Histogram of Precipitation (prcp)", xlab="prcp")
hist(df_advanced2$tmp, main="Histogram of Temperature (tmp)", xlab="tmp")

# # Check collinearity using the VIF
# vif_check <- performance::check_collinearity(df_advanced2)
# print(vif_check)
# 


library(corrplot)
numeric_vars <- df_advanced2[, sapply(df_advanced2, is.numeric)]
cor_matrix <- cor(na.omit(numeric_vars))
corrplot(cor_matrix)

ggplot(df_advanced2, aes(x = richness)) + geom_histogram(bins = 20) + facet_wrap(~ advanced)

# all predictors, no interactions
adv_basic_gam <- gam(advanced ~ s(prcp, k = 15) + s(tmp, k = 15) + s(spei12, k = 15) + 
                   s(distance_edge, k = 15) +
                   s(depth_extract, k =15) +
                   s(disturbance_severity) +
                   s(clay_extract) +
                   s(av.nitro, k =15) +
                   s(management_intensity, k = 10) +
                   s(country_pooled, bs = 're') # + clim_grid + clim_class
                 , 
                 family = binomial(link = "logit"), data = df_advanced2, method = "REML")
summary(adv_basic_gam)
gam.check(adv_basic_gam)
appraise(adv_basic_gam)
windows()
plot.gam(adv_basic_gam, page = 1)


# Model 1: Add interaction between prcp and tmp
adv_interaction_model_1 <- gam(advanced ~ s(prcp, k = 15) + s(tmp, k = 15) + s(spei12, k = 15) + 
                             s(distance_edge) + s(depth_extract) + s(disturbance_severity) + 
                             s(clay_extract) + s(av.nitro) +
                               s(management_intensity, by = country_pooled, bs = "re") +
                             s(country_pooled, bs = "re") +
                             ti(prcp, tmp, k = 15) +# Add interaction
                               s(x,y) +
                               s(clim_grid, bs = "re"), 
                           family = binomial(link = "logit"), 
                           data = df_advanced2, 
                           method = "REML")

# Model 2: Add interaction between prcp and spei12
adv_interaction_model_2 <- gam(advanced ~ s(prcp, k = 15) + s(tmp, k = 15) + s(spei12, k = 15) + 
                             s(distance_edge) + s(depth_extract) + s(disturbance_severity) + 
                             s(clay_extract) + s(av.nitro) +
                               #s(management_intensity, k = 5) + 
                               s(management_intensity, by = country_pooled, bs = "re") +
                               s(country_pooled, bs = "re") +
                             ti(prcp, spei12, k = 5), # Add interaction
                           family = binomial(link = "logit"), 
                           data = df_advanced2, 
                           method = "REML")

# Model 3: Add interaction between tmp and spei12
adv_interaction_model_3 <- gam(advanced ~ s(prcp, k = 15) + s(tmp, k = 15) + s(spei12, k = 15) + 
                             s(distance_edge) + s(depth_extract) + s(disturbance_severity) + 
                             s(clay_extract) + s(av.nitro) + 
                               s(management_intensity, by = country_pooled, bs = "re") +
                               #s(management_intensity, k = 5) + 
                             s(country_pooled, bs = "re") +
                             ti(tmp, spei12, k = 5) , # Add interaction
                           family = binomial(link = "logit"), 
                           data = df_advanced2, 
                           method = "REML")


# include the management vary by country (random slopes), add clim_grid and s(x,y)
adv_interaction_model_3_upd <- gam(advanced ~ s(prcp, k = 15) + 
                                     s(tmp, k = 15) + 
                                     s(spei12, k = 15) + 
                                   s(distance_edge) + 
                                   s(depth_extract) + 
                                     s(disturbance_severity) + 
                                 s(clay_extract) + 
                                   s(av.nitro) +
                                   s(management_intensity, by = country_pooled, bs = "re") +
                                   #s(management_intensity, k = 5) + 
                                   s(clim_grid, bs = "re") +
                                   s(x, y) +
                                 s(country_pooled, bs = "re") +
                                 ti(tmp, spei12, k = 5) , # Add interaction
                               family = binomial(link = "logit"), 
                               data = df_advanced2, 
                               method = "REML")

adv_interaction_model_3_upd2 <- gam(advanced ~ s(prcp, k = 15) + 
                                     s(tmp, k = 15) + 
                                     s(spei12, k = 15) + 
                                     s(distance_edge) + 
                                     s(depth_extract) + 
                                     s(disturbance_severity) + 
                                     s(clay_extract) + 
                                     s(av.nitro) +
                                     s(management_intensity, by = country_pooled, bs = "re") +
                                     #s(management_intensity, k = 5) + 
                                     s(clim_grid, bs = "re") +
                                     s(x, y) +
                                     s(country_pooled, bs = "re") +
                                      ti(tmp, spei12, k = 5) + #, # Add interaction
                                      ti(tmp, prcp, k = 5),
                                   family = binomial(link = "logit"), 
                                   data = df_advanced2, 
                                   method = "REML")


adv_interaction_model_3_upd3 <- gam(advanced ~ s(prcp_c, k = 15) + 
                                      s(tmp_c, k = 15) + 
                                      s(spei12, k = 15) + 
                                      s(distance_edge, k = 20) + 
                                      s(depth_extract) + 
                                      s(disturbance_severity) + 
                                      s(clay_extract) + 
                                      s(av.nitro) +
                                      s(management_intensity, by = country_pooled, bs = "re", k = 20) +
                                      #s(management_intensity, k = 5) + 
                                      s(clim_grid, bs = "re") +
                                      s(x, y) +
                                      s(country_pooled, bs = "re") +
                                      #ti(management_intensity, disturbance_severity) +
                                      #ti(tmp_c, spei12, k = 15) + #, # Add interaction
                                      ti(tmp_c, prcp_c, k = 15),
                                    family = binomial(link = "logit"), 
                                    data = df_advanced2, 
                                    method = "REML")


adv_interaction_model_3_upd4 <- gam(advanced ~ s(prcp, k = 15) + 
                                      s(tmp, k = 15) + 
                                      s(spei12, k = 15) + 
                                      s(distance_edge, k = 20) + 
                                      s(depth_extract) + 
                                      s(disturbance_severity) + 
                                      s(clay_extract) + 
                                      s(av.nitro) +
                                      s(management_intensity, by = country_pooled, bs = "re", k = 20) +
                                      #s(management_intensity, k = 5) + 
                                      s(clim_grid, bs = "re") +
                                      s(x, y) +
                                      s(country_pooled, bs = "re")# +
                                      #ti(management_intensity, disturbance_severity) +
                                      #ti(tmp, spei12, k = 5) + #, # Add interaction
                                      #ti(tmp, prcp, k = 5)
                                      ,
                                    family = binomial(link = "logit"), 
                                    data = df_advanced2, 
                                    method = "REML")


# remove clim_grid
adv_interaction_model_3_upd5 <- gam(advanced ~ s(prcp, k = 15) + 
                                      s(tmp, k = 15) + 
                                      s(spei12, k = 15) + 
                                      s(distance_edge, k = 20) + 
                                      s(depth_extract) + 
                                      s(disturbance_severity) + 
                                      s(clay_extract) + 
                                      s(av.nitro) +
                                      s(management_intensity, by = country_pooled, bs = "re", k = 20) +
                                      #s(management_intensity, k = 5) + 
                                      #s(clim_grid, bs = "re") +
                                      s(x, y) +
                                      s(country_pooled, bs = "re")# +
                                    #ti(management_intensity, disturbance_severity) +
                                    #ti(tmp, spei12, k = 5) + #, # Add interaction
                                    #ti(tmp, prcp, k = 5)
                                    ,
                                    family = binomial(link = "logit"), 
                                    data = df_advanced2, 
                                    method = "REML")




AIC(adv_interaction_model_3_upd5,
  adv_interaction_model_3_upd4, 
  adv_interaction_model_3_upd3,
     adv_interaction_model_3_upd2,
     adv_interaction_model_3_upd, 
     adv_interaction_model_3,
     adv_interaction_model_2,
     adv_interaction_model_1)

BIC( adv_interaction_model_3_upd5,
  adv_interaction_model_3_upd4,
  adv_interaction_model_3_upd3,
     adv_interaction_model_3_upd2,
     adv_interaction_model_3_upd, 
     adv_interaction_model_3,
     adv_interaction_model_2,
     adv_interaction_model_1)


predicted_interaction <- ggpredict(adv_interaction_model_3_upd3, terms = c("tmp_c", "prcp_c"))
plot(predicted_interaction)

# ANOVA to see if tmp differs by clim_grid
anova_tmp <- aov(tmp ~ clim_grid, data = df_advanced2)
summary(anova_tmp)

# ANOVA to see if prcp differs by clim_grid
anova_prcp <- aov(prcp ~ clim_grid, data = df_advanced2)
summary(anova_prcp)

# Boxplot for tmp across clim_grid levels
ggplot(df_advanced2, aes(x = clim_grid, y = tmp)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Temperature Variation Across Clim Grid Levels", x = "Clim Grid", y = "Temperature (tmp)")

# Boxplot for prcp across clim_grid levels
ggplot(df_advanced2, aes(x = clim_grid, y = prcp)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Precipitation Variation Across Clim Grid Levels", x = "Clim Grid", y = "Precipitation ")

windows()
plot.gam(adv_interaction_model_3_upd2, page = 1)
appraise(adv_interaction_model_3_upd2)
summary(adv_interaction_model_3_upd2)
gam.check(adv_interaction_model_3_upd2)
k.check(adv_interaction_model_3_upd2)

hist(df_advanced2$distance_edge)

# see influential points:
influence_data <- influence.gam(adv_interaction_model_3_upd2)
summary(influence_data)

# identify observations with high influence
high_influence <- which(influence_data > 0.1)  # Adjust the threshold as needed
df_advanced2[high_influence, ]


draw(adv_interaction_model_3_upd2)

fin.m.advanced <- adv_interaction_model_3_upd2# adv_interaction_model_3_upd3
y_lab = "Probability Advanced"

# Test for spatial autocorrelation
model_residuals <- residuals(fin.m.advanced, type = "pearson")

# Load necessary libraries
library(spdep)

# Create coordinates matrix
coords <- cbind(df_advanced2$x, df_advanced2$y)

# Create a spatial neighbors object (e.g., using k-nearest neighbors)
# Adjust k based on the density and distribution of your data points
nb <- knn2nb(knearneigh(coords, k = 5))

# Convert neighbors list to a weights list
listw <- nb2listw(nb, style = "W")
# Perform Moran's I test
moran_test <- moran.test(model_residuals, listw)
print(moran_test)
# 0.57



#### Plots advanced----------------------------
summary(fin.m.advanced)
y_lab = 'Prob adv'

# Predict management intensity effects by country
predicted_mgmt_intensity <- ggpredict(fin.m.advanced, terms = c("management_intensity", "country_pooled"))
# Generate predictions for stem regeneration across different countries
predicted_countries <- ggpredict(fin.m.advanced, terms = "country_pooled")

# Create plots for individual variables
plot_prcp   <- create_plot(fin.m.advanced, "prcp", df_advanced2 , "Prcp 0.6", scatter_y = "advanced",line_color = "blue", fill_color = "blue")
plot_tmp    <- create_plot(fin.m.advanced, "tmp", df_advanced2, "Tmp 0.9", scatter_y = "advanced",line_color = "red", fill_color = "red")
plot_spei12 <- create_plot(fin.m.advanced, "spei12", df_advanced2, "spei12 0.07", scatter_y = "advanced",line_color = "darkgreen", fill_color = "green")
plot_disturbance_severity <- create_plot(fin.m.advanced, "disturbance_severity", df_advanced2, "Dist Sev 0.4", scatter_y = "advanced", line_color = "purple", fill_color = "purple")
plot_distance_edge <- create_plot(fin.m.advanced, "distance_edge", df_advanced2, "Dist edge 0.3", x_limit = c(0, 600), scatter_y = "advanced", line_color = "grey", fill_color = "grey")
#plot_management_intensity <- create_plot(fin.m.advanced, "management_intensity", df_advanced2, "Mng intend",scatter_y = "delayed", line_color = "grey", fill_color = "grey")


# Create interaction plots
plot_interaction1 <- create_interaction_plot(fin.m.advanced, 
                                             c("prcp", "tmp"), "prcp & tmp *", df_advanced2) +
  guides(color = guide_legend(title = "tmp"), fill = guide_legend(title = "tmp")) + # Rename legend to 'tmp'
  theme(legend.position.inside = c(0.7, 0.7),
        #legend.position.inside = "top right", # Place the legend in the upper right corner inside the plot
        legend.background = element_rect(fill = "white", color = NA), # Add a white background
        legend.title = element_text(size = 10), # Adjust the title size
        legend.text = element_text(size = 8))   # Adjust the text size

plot_interaction2 <- create_interaction_plot(fin.m.advanced, c("prcp", "spei12"), "precp & spei 0.17", 
                                             df_advanced2) +
  guides(color = guide_legend(title = "tmp"), fill = guide_legend(title = "tmp")) + # Rename legend to 'tmp'
  theme(legend.position.inside =  c(0.7, 0.7), # Place the legend in the upper right corner inside the plot
        legend.background = element_rect(fill = "white", color = NA), # Add a white background
        legend.title = element_text(size = 10), # Adjust the title size
        legend.text = element_text(size = 8))   # Adjust the text size

# Plot predicted effects by country
plot_countries <- ggplot(predicted_countries, aes(x = x, y = predicted, fill = x)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = "", x = "Country", y = y_lab) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the effect of management intensity on stem regeneration by country
plot_mgmt_intensity <- ggplot(predicted_mgmt_intensity, aes(x = x, y = predicted, group = group)) +
  geom_line(aes(color = group), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  labs(title = "", 
       x = "Management Intensity", y = y_lab, color = "Country", fill = "Country") +
  theme_minimal() +
  facet_wrap(~group) +
  theme(legend.position = "none")


# !!!
library(cowplot)
#plot_grid(p1, p2, labels = "auto", label_size = 12)
comb_upper  <- plot_grid(plot_interaction1, plot_interaction2,ncol = 2)
comb_lower <- plot_grid(plot_countries,plot_prcp, plot_tmp, plot_spei12, plot_disturbance_severity,
                          plot_distance_edge, ncol = 2)
# Combine the individual plots into one figure
combined_plot <- plot_grid(comb_upper,
                           comb_lower,
                           ncol = 1, nrow = 2,
                           rel_heights = c(1, 3))

combined_plot


ggsave('outFigs/fig_advanced_drivers.png', plot = combined_plot, width =5,
       height = 7,bg = 'white')



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
                  dv.labels = paste0("Explained Deviance: ", round(100 * summary(fin.m.advanced)$dev.expl, 2), "%"), 
                  file = "outTable/full_drivers_advanced.doc")











## Drivers regeneration density pooled ---------------------------------------------


# Check the structure of the data
str(df_stem_regeneration2)

# Check for missing values
sum(is.na(df_stem_regeneration2))

# Inspect the distribution of each predictor
hist(df_stem_regeneration2$prcp, main="Histogram of Precipitation (prcp)", xlab="prcp")
hist(df_stem_regeneration2$tmp, main="Histogram of Temperature (tmp)", xlab="tmp")

# Check collinearity using the VIF
vif_check <- performance::check_collinearity(df_stem_regeneration2)
print(vif_check)



library(corrplot)
numeric_vars <- df_stem_regeneration2[, sapply(df_stem_regeneration2, is.numeric)]
cor_matrix <- cor(na.omit(numeric_vars))
corrplot(cor_matrix)

ggplot(df_stem_regeneration2, aes(x = richness)) + geom_histogram(bins = 20) + facet_wrap(~ stem_regeneration    )

basic_gam <- gam(stem_regeneration     ~ s(prcp, k = 15) + s(tmp, k = 15) + s(spei12, k = 15) + 
                   s(distance_edge) +
                   s(depth_extract, k = 15) +
                   s(disturbance_severity) +
                   s(clay_extract) +
                   s(av.nitro, k =20) +
                   s(management_intensity, k = 10) +
                   s(country_pooled, bs = 're') # + clim_grid + clim_class
                 , 
                 family = Tweedie(p=1.5), data = df_stem_regeneration2, method = "REML")
summary(basic_gam)
gam.check(basic_gam)
appraise(basic_gam)
windows()
plot.gam(basic_gam, page = 1)


# Model 1: Adding interaction between prcp and tmp
interaction_model_1 <- gam(stem_regeneration ~ s(prcp, k = 15) + s(tmp, k = 15) + 
                             s(spei12, k = 15) + s(distance_edge) + s(depth_extract, k = 15) + 
                             s(disturbance_severity) + s(clay_extract) + 
                             s(av.nitro, k = 20) + s(management_intensity, k = 10) + 
                             s(country_pooled, bs = "re") + ti(prcp, tmp, k = 10),
                           family = tw(), method = "REML", data = df_stem_regeneration2)

summary(interaction_model_1)

# Model 2: Adding interaction between spei12 and tmp
interaction_model_2 <- gam(stem_regeneration ~ s(prcp, k = 15) + s(tmp, k = 15) + 
                             s(spei12, k = 15) + s(distance_edge) + s(depth_extract, k = 15) + 
                             s(disturbance_severity) + s(clay_extract) + 
                             s(av.nitro, k = 20) + s(management_intensity, k = 10) + 
                             s(country_pooled, bs = "re") + ti(spei12, tmp, k = 10),
                           family = tw(), method = "REML", data = df_stem_regeneration2)

summary(interaction_model_2)

# Model 3: Adding interaction between prcp and management_intensity
interaction_model_3 <- gam(stem_regeneration ~ s(prcp, k = 15) + s(tmp, k = 15) + 
                             s(spei12, k = 15) + s(distance_edge) + s(depth_extract, k = 15) + 
                             s(disturbance_severity) + s(clay_extract) + 
                             s(av.nitro, k = 20) + s(management_intensity, k = 10) + 
                             s(country_pooled, bs = "re") + ti(prcp, management_intensity, k = 10),
                           family = tw(), method = "REML", data = df_stem_regeneration2)

summary(interaction_model_3)

# Model 4: Multiple interactions & account for clim_grid to avoid spatial autocorrelation
interaction_model_4 <- gam(stem_regeneration ~ s(prcp, k = 15) + s(tmp, k = 15) + 
                             s(spei12, k = 15) + s(distance_edge) + s(depth_extract, k = 15) + 
                             s(disturbance_severity) + s(clay_extract) + 
                             s(av.nitro, k = 20) + s(management_intensity,by = country_pooled, k = 10) + 
                             s(country_pooled, bs = "re") + ti(prcp, tmp, k = 10) + 
                             ti(spei12, tmp, k = 10) +
                             s(x,y) +
                             s(clim_grid, bs = "re") +
                             country_pooled
                             #ti(prcp, management_intensity, k = 10)
                             ,
                           family = tw(), method = "REML", data = df_stem_regeneration2)

summary(interaction_model_4)


# Modify the model to use a random slope for management_intensity by country
interaction_model_5 <- gam(
  stem_regeneration ~ s(prcp, k = 15) + s(tmp, k = 15) + s(spei12, k = 15) +
    s(distance_edge) + s(depth_extract, k = 15) + s(disturbance_severity) +
    s(clay_extract) + s(av.nitro, k = 20) +
    # Change from separate smooths to a random slope model for management_intensity by country
    s(management_intensity, by = country_pooled, bs = "re") +
    s(country_pooled, bs = "re") +  # Random intercept for country
    ti(prcp, tmp, k = 10) + ti(spei12, tmp, k = 10) +
    s(x, y) + s(clim_grid, bs = "re"), 
  family = tw(), method = "REML", data = df_stem_regeneration2
)


# View the summary of the updated model
summary(interaction_model_5)
summary(interaction_model_4)
appraise(interaction_model_5)

# store the best model for regeneration density
fin.m.reg.density <- interaction_model_5 

# test for autocorrelation

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


# Compare models
AIC(basic_gam, interaction_model_1, interaction_model_2, interaction_model_3, interaction_model_4,
    interaction_model_5,
    interaction_model_random_smooth,interaction_model_updated)

# Display model summaries to understand which terms are significant
summary(basic_gam)
summary(interaction_model_1)
summary(interaction_model_2)
summary(interaction_model_3)
summary(interaction_model_4)
summary(interaction_model_5)


# Update model to include interaction between management_intensity and country_pooled
interaction_model_updated <- gam(stem_regeneration ~ 
                                   s(prcp, k = 15) + s(tmp, k = 15) + s(spei12, k = 15) + 
                                   s(distance_edge) + s(depth_extract, k = 15) + s(disturbance_severity) + 
                                   s(clay_extract) + s(av.nitro, k = 20) + 
                                   ti(management_intensity, country_pooled, bs = "fs") + 
                                   ti(prcp, tmp, k = 10) + ti(spei12, tmp, k = 10) + 
                                   s(x,y) +
                                   country_pooled, 
                                 family = tw(), method = "REML", data = df_stem_regeneration2)
summary(interaction_model_updated)


# Adding random smooth effect for country
interaction_model_random_smooth <- gam(stem_regeneration ~ 
                                         s(prcp, k = 15) + s(tmp, k = 15) + s(spei12, k = 15) + 
                                         s(distance_edge) + s(depth_extract, k = 15) + s(disturbance_severity) + 
                                         s(clay_extract) + s(av.nitro, k = 20) + 
                                         s(management_intensity, by = country_pooled) + 
                                         ti(prcp, tmp, k = 10) + ti(spei12, tmp, k = 10) + 
                                         s(x,y) +
                                         s(country_pooled, bs = "re"),
                                       family = tw(), method = "REML", data = df_stem_regeneration2)
summary(interaction_model_random_smooth)





# PLot 
y_lab = 'Stem dens [#/ha]'

summary(fin.m.reg.density)
# Predict management intensity effects by country
predicted_mgmt_intensity <- ggpredict(fin.m.reg.density, terms = c("management_intensity", "country_pooled"))
# Generate predictions for stem regeneration across different countries
predicted_countries <- ggpredict(fin.m.reg.density, terms = "country_pooled")


# Create plots for individual variables
plot_prcp <- create_plot(fin.m.reg.density, "prcp", df_stem_regeneration2, "prcp**", line_color = "blue", fill_color = "blue")
plot_tmp <- create_plot(fin.m.reg.density, "tmp", df_stem_regeneration2, "tmp**", line_color = "red", fill_color = "red")
plot_spei12 <- create_plot(fin.m.reg.density, "spei12", df_stem_regeneration2, "spei12 0.9", line_color = "darkgreen", fill_color = "green")
plot_disturbance_severity <- create_plot(fin.m.reg.density, "disturbance_severity", df_stem_regeneration2, "Dist Severity*", line_color = "purple", fill_color = "purple")
plot_distance_edge <- create_plot(fin.m.reg.density, "distance_edge", df_stem_regeneration2, "Dist edge 0.12", x_limit = c(0, 600), line_color = "grey", fill_color = "grey")
plot_management_intensity <- create_plot(fin.m.reg.density, "management_intensity", df_stem_regeneration2, "Management_intensity", line_color = "grey", fill_color = "grey")


# Create interaction plots
plot_interaction1 <- create_interaction_plot(fin.m.reg.density, 
                                             c("prcp", "tmp[8,9,10]"), "prcp & tmp 0.05", df_stem_regeneration2) +
  labs(color = "tmp", fill = "tmp") 
plot_interaction2 <- create_interaction_plot(fin.m.reg.density, 
                                             c("spei12", "tmp[8,9,10]"), "spei12 & tmp 0.5",
                                             df_stem_regeneration2) +labs(color = "tmp", fill = "tmp") 

# Plot predicted effects by country
plot_countries <- ggplot(predicted_countries, aes(x = x, y = predicted, fill = x)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(width = 0.9)) +
  labs(title = "", x = "Country", y = y_lab) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the effect of management intensity on stem regeneration by country
plot_mgmt_intensity <- ggplot(predicted_mgmt_intensity, aes(x = x, y = predicted, group = group)) +
  geom_line(aes(color = group), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  labs(title = "", 
       x = "Manag Intensity", y = y_lab, color = "Country", fill = "Country") +
  theme_minimal() +
  facet_wrap(~group) +
  theme(legend.position = "none")


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











## Wilcox PLOTS: sites differences: delayed vs advaced ----------------------------

# Select relevant columns including 'RegenerationStatus' and the desired variables
variables_to_plot <- c(
  "rIVI", "sum_stems_mature", "prcp", 
  "richness",
  "n_vertical",
 # "spei1", "drought_spei1",                      
  "spei3", #"drought_spei3",
                        "spei12", #"drought_spei12", 
                      
                      # "spei24", "drought_spei24" ,
                       "tmp", "av.nitro", "depth_extract", "management_intensity", 
 "distance_edge" ,
 "salvage_intensity",
 "clay_extract",
 "clay_extract",
 "disturbance_severity",
                       "adv_delayed"
                       )

#                       "country_pooled", "clim_grid", "clim_class"

# Step 2: Calculate median and IQR for each variable by regeneration status
summary_stats <- 
  df_fin %>%
  na.omit() %>% 
    dplyr::select(all_of(variables_to_plot)) %>% 
  gather(key = "Variable", value = "Value", -adv_delayed) %>%
  group_by(adv_delayed, Variable) %>%
  summarise(
    Median = median(Value, na.rm = TRUE),
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE)
  ) %>%
  mutate(IQR_Lower = Median - Q1, 
         IQR_Upper = Q3 - Median)

# Step 3: Create the median and IQR plot using ggplot2
ggplot(summary_stats, aes(x = Variable, y = Median, color = adv_delayed)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Points for median
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2, 
                position = position_dodge(width = 0.5)) +          # Error bars for IQR
  labs(title = "Median and IQR Plot for Delayed vs Advanced Regeneration",
       x = "Variables", y = "Median and IQR") +
  theme_minimal() +
  facet_wrap(.~Variable, scales = 'free')+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("blue", "red")) +  # Set colors for delayed vs advanced
  theme(legend.title = element_blank())            # Remove legend title


# Step 1: Reshape the data to long format
df_long <- df_fin %>%
  na.omit() %>% 
  dplyr::select(all_of(variables_to_plot), adv_delayed) %>% 
  gather(key = "Variable", value = "Value", -adv_delayed)

# Calculate max value for each variable to set label.y properly
df_long_summary <- df_long %>%
  group_by(Variable) %>%
  summarize(max_value = max(Value, na.rm = TRUE)) %>%
  mutate(label_y = max_value + (0.1 * max_value)) # Add a 10% margin above the max value for label

# Create a named vector for easy lookup in stat_compare_means
label_y_positions <- setNames(df_long_summary$label_y, df_long_summary$Variable)

# Plot using ggboxplot
# Plot using ggboxplot
p_boxplot_wilcox <- ggboxplot(df_long, x = "adv_delayed", y = "Value", 
                              color = "adv_delayed", palette = c("blue", "red"),
                              facet.by = "Variable", scales = "free_y", 
                              ylab = "Values", xlab = "Regeneration Status") +
  stat_compare_means(aes(group = adv_delayed), method = "wilcox.test", 
                     label = "p.signif", 
                     label.y = label_y_positions[as.character(df_long$Variable)], # Use the calculated y positions
                     size = 4,
                     label.x = 1.5) +  # Position labels between the two boxplots
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(title = "",
       x = "Reg. Status", y = "Vals")

p_boxplot_wilcox

# Save the combined plot (optional)
ggsave("outFigs/p_boxplot_wilcox.png", p_boxplot_wilcox, width = 6, height = 8, 
       dpi = 300,  bg = 'white')




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



#### Stem density vs IVI -----------------------------------------------------


# Save models --------------------------------------------------

# Save the model object and input data
save(#fin.m.delayed, df_delayed2, 
     fin.m.advanced, df_advanced2, 
     #fin.m.advanced_juv, df_advanced_juv2,
     #fin.m.reg.density, df_stem_regeneration2,
     file = "outData/stem_density_models.RData")

