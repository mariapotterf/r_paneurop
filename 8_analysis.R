# Analyse the data based on final table: with several speis

# run cluster analysis for climate& environmnet
# for structural data
# get summary tables


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

source('my_functions.R')

# Read data -----------------------------------------------------------------------

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
  left_join(xy_site, by = join_by(site))


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
          clim_class = as.factor(clim_class))


# export as xy coordinates climate cluster categorized 
df_fin_clim_clust_xy <- st_as_sf(df_fin, coords = c("x", "y"), crs = crs(xy))  # Replace '4326' with the appropriate CRS if known

# Step 2: Check the structure of the sf object to ensure everything is correct
print(df_fin_clim_clust_xy)

# Step 3: Export the data as a GeoPackage (GPKG)
st_write(df_fin_clim_clust_xy, "outData/xy_clim_cluster.gpkg", layer = "df_fin", driver = "GPKG")




# Analysis ------------------------------------------------------------
# desription of current regeneration state
# how many clusters do I have by managemnet types?
# how do vegetation variables behave along management cluster?

# What is the current state of teh forest regeneration?
#  stucture? 
#    - stem density
#    - vertical structure
#  Composition? 
#    - richness
#    - species dominance
#  along managemnet gradient? 
# H0: presence of regeneration, dominance of teh avanced regeneration
# how often classes cooccur? presence/absence of vertical clases
# higher management intensity, higher species richness
# higher managemnet, less of teh advanced regeneration - 
#   due to the salvage logging and following 

# predictors: investigate, which SPEI is the best? 
#             correlate with stem density

# drivers:
# - apply drivers for all 4 veg indicators??
# effect of climate: 
#   - higher drought, less regeneration; higher change of delayed regeneration?

# effect of disturbance size:
#   - higher patch size, less regeneration
#   - higher severity, low regenerations

# environmental condistins;
#  soil: 
#    - more stem density at better soil conditions (higher depth, nutrient content)






# Make 2d denisty plots: 

# Merge structure & composition into single space 


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
variables <- c("tmp", "prec", "tmp_z", "prcp_z", "spei1", "spei3", "spei6", "spei12", "spei24")

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



# try simple glm: values between 0-1 - use beta distribution:
# # make sure that no dat are equal 0 or 1
# df_fin
# # Ensure that the transformed variable is strictly between 0 and 1, not 0 or 1
# dat_fin$tr_agg_doy  <- pmin(pmax(dat_fin$tr_agg_doy, 1e-4),  1 - 1e-4)
# dat_fin$tr_peak_doy <- pmin(pmax(dat_fin$tr_peak_doy, 1e-4), 1 - 1e-4)
# 
# 
# model_lag0 <- glmmTMB(model_lag0_formula, family = beta_family(link = "logit"), data = data, 
#                       na.action = na.exclude)
# 
# # Assuming your model is an 'lme4' or 'glmmTMB' object
# (vif_values <- performance::check_collinearity(m1))




###### Climate space:  density plot with raster -------------------------------

library(MASS)
library(RColorBrewer)

blue_colors <- brewer.pal(5, "Greens")

# Add 'white' at the beginning
my_colors <- c("white", blue_colors)

# Function to calculate levels
calculate_levels <- function(density, percentages) {
  sorted_density <- sort(density, decreasing = TRUE)
  cum_density <- cumsum(sorted_density) / sum(sorted_density)
  sapply(percentages, function(p) {
    sorted_density[max(which(cum_density <= p))]
  })
}



##### Define the function to prepare density data --------------------------
prepare_density_data <- function(df_fin, 
                                 x = 'rIVI', y = 'stem_density',
                                 percentages = c(0.25, 0.50, 0.75, 0.90, 1)) {
  # Calculate 2D density
  d <- kde2d(df_fin[[x]], df_fin[[y]], n = 300)
  density_values <- d$z
  
  # Set densities outside the range of your data to zero
  density_values[density_values < 1e-6] <- 0
  
  # Calculate the levels for specified percentages
  levels <- calculate_levels(as.vector(density_values), percentages)
  
  # Prepare data for ggplot
  plot_data <- expand.grid(x = d$x, y = d$y) %>%
    mutate(density = as.vector(density_values))
  
  # Use cut to create factor levels, including one for zero density
  plot_data$level <- cut(plot_data$density, breaks = c(-Inf, levels, Inf), labels = FALSE, include.lowest = TRUE)
  
  return(plot_data)
}

# make function to produce  2D desity plot
make_2D_plot <- function(data = plot_data, 
                         x_var = 'x', # Default column names, can be changed
                         y_var = 'y', # Default column names, can be changed
                         x_lab = "X Axis", # Default axis labels, can be customized
                         y_lab = "Y Axis") {#, # Default axis labels, can be customized) 
  
  # Convert the y_var argument to a symbol to use in aes()
  x_var_sym <- rlang::sym(x_var)
  y_var_sym <- rlang::sym(y_var)
  
  p <- ggplot(data) +
    geom_raster(aes(x = !!x_var_sym, 
                    y = !!y_var_sym, 
                    fill = factor(level)), 
                alpha = 0.8) +
    scale_fill_manual(values = my_colors,
                      labels = c("", rev(c("25%", "50%", "75%", "90%", "100%"))),
                      name = "") +
    ylab(y_lab) +
    xlab(x_lab) +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          aspect.ratio = 1,
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.3, "cm"),
          legend.background = element_blank()) #+
  
  return(p)
}

##### Density plots: each variable against another one --------------------
# 1. rIVI vs stem_density
# 2. rIVI vs richness
# 3. rIVI vs n_vertical
# 4. stem_density vs richness
# 5. stem_density vs n_vertical
# 6. richness vs n_vertical 


######### rIVI vs sum_stems 
df_plot_data <- prepare_density_data(df_fin, 
                                     x = 'stem_density',
                                     y = 'rIVI')


p_IVI_stems <- make_2D_plot(df_plot_data, 
                            x_lab = "Stem density [n/ha]",
                            y_lab = "Importance value [%]")

(p_IVI_stems)
####### rIVI vs richnes
df_plot_data <- prepare_density_data(df_fin, 
                                     x = 'rIVI', 
                                     y = 'richness')


p_IVI_richness <- make_2D_plot(df_plot_data, 
                               x_lab = "Importance value [%]", # Default axis labels, can be customized
                               y_lab = "Richness [counts]"  )


####### rIVI vs vertical 
df_plot_data <- prepare_density_data(df_fin, 
                                     x = 'rIVI', 
                                     y = 'n_vertical')


p_IVI_vert <- make_2D_plot(df_plot_data, 
                           x_lab = "Importance value [%]", # Default axis labels, can be customized
                           y_lab = "# vertical layers [counts]"  )



####### stem_density vs richness 
df_plot_data <- prepare_density_data(df_fin, 
                                     x = 'stem_density', 
                                     y = 'richness')


p_dens_richness <- make_2D_plot(df_plot_data, 
                                x_lab = "Stem density [n/ha]", # Default axis labels, can be customized
                                y_lab = "Richness [counts]"  )



####### stem_density vs n_vertcal 
df_plot_data <- prepare_density_data(df_fin, 
                                     x = 'stem_density', 
                                     y = 'n_vertical')


p_dens_vertical <- make_2D_plot(df_plot_data, 
                                x_lab = "Stem density [n/ha]", # Default axis labels, can be customized
                                y_lab = "# vertical layers [counts]"  )


####### n_vertical vs richness 
df_plot_data <- prepare_density_data(df_fin, 
                                     x = 'n_vertical', 
                                     y = 'richness')


p_vertical_richness <- make_2D_plot(df_plot_data, 
                                    x_lab = "# vertical layers [counts]", # Default axis labels, can be customized
                                    y_lab = "Richness [counts]"  )



ggarrange(p_IVI_stems,
          p_IVI_vert,
          p_IVI_richness,
          p_dens_richness,
          p_dens_vertical,
          p_vertical_richness, nrow = 2, ncol = 3, common.legend = T, legend = 'bottom', align = 'hv',
          labels = c("[a]","[b]","[c]","[d]","[e]","[f]"))
#### DEnsity plot wth management 
###### rIVI vs management 
df_plot_data <- prepare_density_data(df_fin, 
                                     y = 'rIVI', 
                                     x = 'management_intensity')


p_IVI_manag <- make_2D_plot(df_plot_data, 
                            y_lab = "Importance value [%]", # Default axis labels, can be customized
                            x_lab = "Management intensity [%]"  )

(p_IVI_manag)

###### richness vs management 
df_plot_data <- prepare_density_data(df_fin, 
                                     y = 'richness', 
                                     x = 'management_intensity')


p_richness_manag <- make_2D_plot(df_plot_data, 
                                 y_lab = "Richness [#]", # Default axis labels, can be customized
                                 x_lab = "Management intensity [%]"  )

(p_richness_manag)



###### stem_density vs management
df_plot_data <- prepare_density_data(df_fin, 
                                     y = 'stem_density', 
                                     x = 'management_intensity')


p_stem_dens_manag <- make_2D_plot(df_plot_data, 
                                  y_lab = "Stem density [n/ha]", # Default axis labels, can be customized
                                  x_lab = "Management intensity [%]"  )

(p_stem_dens_manag)




###### n_vertical vs management 
df_plot_data <- prepare_density_data(df_fin, 
                                     y = 'n_vertical', 
                                     x = 'management_intensity')


p_vertical_manag <- make_2D_plot(df_plot_data, 
                                 y_lab = "Vertical layers [#]", # Default axis labels, can be customized
                                 x_lab = "Management intensity [%]"  )

(p_vertical_manag)



# Merge plots indicators vs management intensity 

ggarrange(p_IVI_manag, p_richness_manag, 
          p_stem_dens_manag, p_vertical_manag, 
          nrow = 2, ncol = 2,
          common.legend = T, #legend.position = 'right', 
          align = 'hv',
          labels = c("[a]", "[b]","[c]","[d]"))




# 2. species composition per stems: saplings vs juveniles?? --------------
# how many plots with Mature trees?
# 

# get overall number for the species compositions:
stem_dens_species_long %>% 
  group_by(Species) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup(.) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share       = sum_stems/sum_vegType*100) %>% 
  arrange(desc(share)) #%>% 
#View()

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




dom_species <- stem_dens_species_long %>% 
  filter(VegType != 'Mature' ) %>%  # Matures: on 191
  group_by(Species, VegType) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup(.) %>% 
  group_by(VegType) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share       = sum_stems/sum_vegType*100) #%>% # 'share by species' is calculated within the vertical type group 

table(dom_species$Species,dom_species$VegType)

# dominant species by country?
dom_species_country <- 
  stem_dens_species_long %>% 
  filter(VegType != 'Mature' ) %>%  # Matures: on 191
  group_by(Species, VegType, country ) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup(.) %>% 
  group_by(VegType,country) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share       = sum_stems/sum_vegType*100)



# SElect the 5 species with teh highest share per management, vertical structure
top_species_per_country <- dom_species_country %>%
  ungroup(.) %>% 
  group_by(country, VegType) %>%
  arrange(desc(share)) %>%  # Arrange in descending order of share within each group
  slice_max(order_by = share, n = 5)  # Select top 5 species with max share in each group

#View(top_species_per_country)


# SElect the 5 species with teh highest share per vertical structure
top_species_per_group <- dom_species %>%
  group_by(VegType) %>%
  arrange(desc(share)) %>%  # Arrange in descending order of share within each group
  slice_max(order_by = share, n = 10)  # Select top 10 species with max share in each group

(top_species_per_group)

# see all species
top_species_per_group %>% 
  mutate(Species = factor(Species, levels = unique(Species[order(share, decreasing = TRUE)]))) %>%
  # arrange(desc(share)) %>% 
  #filter(
  #dplyr::filter(share > 1) %>% 
  ggplot(aes(x = Species,
             y = share,
             #group = VegType,
             fill = VegType)) +
  geom_bar(stat = "identity", position = position_dodge())+
  labs(xlab = '',
       ylab = 'Stem share/category [%]') +
  facet_grid(.~VegType) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0.5),
        legend.position = 'bottom') #+


#geom_col('')


# Reorder Species within each group based on share
top_species_per_group <- top_species_per_group %>%
  arrange(VegType, desc(share)) %>%
  mutate(Species = factor(Species, levels = unique(Species)),
         VegType = factor(VegType, 
                          levels = c("Saplings", "Juveniles")))


# Create a bar plot with ordered categories
ggplot(top_species_per_group, aes(x = VegType, y = share, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(share), y = share), 
            position = position_dodge(width = 0.9), 
            color = "black", size = 3, vjust = -0.5) + # Adjust vjust for label position
  #  facet_grid(. ~ VegType) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) +
  labs(y = "Share [%]", x = "", fill = "Species") +
  theme_classic()










# Descriptive plots  -----------------------------------------------------------

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





# for climate ------------------------------------------------------------------
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


# test differences between medians of teh classes -----------------------------

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




# For disturbance size and intensity -------------------------------------------

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


# test differences between medians of teh classes -----------------------------

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


# How many vertical layers I have where? --------------------------------------

# Get again individual layers:
df_vert_full <- 
  stem_dens_species_long %>% 
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




#5. test model ---------------------------------


# variables selecteion processs for drivers  ------------------------------------

## predictor selection: SPEI for stem density:  ------------------------------------------------------

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


### check for multicollinearity -----------------------------------------------------

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

### test with univariate models & AIC ----------------------------------------------
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

# Run univariate models for set of dependent variables (stem density)
# and spei predictors to se teh best spei

# List of dependent variables
dependent_vars <- c("sum_stems_juvenile", 
                    "sum_stems_sapling", 
                    "sum_stems_mature"
                    #"stem_density"
)

# List of predictor variables (spei1 to spei24, drought_spei1 to drought_spei24)
predictor_vars <- c("spei1", "spei3", "spei6", "spei12", "spei24",
                    "tmp", "tmp_z", "prcp", "prcp_z", 
                    "drought_spei1", "drought_spei3", "drought_spei6", 
                    "drought_spei12", "drought_spei24",
                    "drought_tmp", "drought_prcp")



# TEST

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
best_predictors <- model_metrics %>%
  mutate(category = case_when(
    grepl("sapl", Dependent  ) ~ "sapling",
    grepl("juven", Dependent  ) ~ "juvenile",
    grepl("mature", Dependent  ) ~ "mature",
    TRUE ~ "other"
  )) %>%
  group_by(Dependent, category) %>%
  slice_min(AIC, n = 3)   # Select the best 3 based on AIC
# slice(which.max(DevianceExplained))

best_predictors

sjPlot::tab_df(model_metrics,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = FALSE,
               file="outTable/find_best_spei.doc",
               digits = 1) 

# we found that for every vegetation vertical cla has different sentiticity to SPEI scale:


### identify the most important predictors -----------------------------------------

# Vector of possible predictors
predictors_all <- c("management_intensity", 
                    # "n_vertical", 
                    #"tmp", 
                    # "prcp", 
                    "tmp_z", 
                    "prcp_z", 
                    # "spei1", 
                    #  "spei3", 
                    #  "spei6", 
                    #  "spei12", 
                    # "spei24", 
                    #"drought_spei1", 
                    "drought_spei3", 
                    #"drought_spei6", 
                    "drought_spei12", 
                    "drought_spei24", 
                    "distance_edge", 
                    "disturbance_severity", 
                    "sand_extract", 
                    "clay_extract", 
                    "depth_extract", 
                    "av.nitro" #, 
                    #"salvage_intensity", 
                    #"protection_intensity" #, 
                    # "clim_cluster_spei3", 
                    #"clim_cluster_spei6" #, 
                    #"region", 
                    #"clim_class", 
                    #"regions", 
                    #"country_abbr", 
                    #"country_pooled"
)



# Get correlation analysis: 

# Calculate correlation matrix for numeric variables
cor_matrix    <- cor(df_fin[, c("sum_stems_juvenile", "sum_stems_sapling", "sum_stems_mature", predictors_all)], method = 'spearman')

#library(corrplot)
# Plot the correlation matrix using corrplot
corrplot(cor_matrix)



# Identify variables using randm forests
library(randomForest)

# Create a Random Forest model to identify important predictors
rf_model <- randomForest(sum_stems_juvenile ~ ., data = df_fin[, c("sum_stems_juvenile", predictors_all)], importance = TRUE)
importance(rf_model)
varImpPlot(rf_model)

# Check correlations among important predictors
important_vars <- df_fin[, c("tmp", "prcp", "tmp_z", "prcp_z", "drought_spei3", "av.nitro")]
cor(important_vars)


# Fit a linear model to check VIF
fit <- lm(sum_stems_juvenile ~ tmp + prcp + tmp_z + prcp_z + drought_spei3 + av.nitro, data = df_fin)

# Check VIF
vif(fit)


# how often does stem density occur?
# 



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
                          sand_extract  +
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
                  sand_extract  +
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
                 disturbance_severity + sand_extract + s(clay_extract, k = 10) +
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

#TEST models ----------------------------------------------------------------------------------





# Use Tweedie family with p=1.5
gam_model <- gam(sum_stems_juvenile ~ s(tmp_z, k = 7) + s(prcp_z, k = 7) + 
                   #s(drought_spei12, k = 10) + 
                   s(distance_edge, k = 10) +
                   disturbance_severity + sand_extract + s(clay_extract, k = 10) +
                   s(depth_extract, k = 10) + s(av.nitro, k = 15) +
                   s(management_intensity, by = country_pooled, k = 10) +  # Country-specific smooth
                   #s(management_intensity , k = 10) +
                   ti(tmp_z, prcp_z) +
                   ti(x,y ) +
                   s(region, bs = 're')+ 
                   s(country_abbr, bs = "re"), 
                 family = Tweedie(p = 1.46),  # Set starting p value
                 data = df_fin)

cor(df_fin$tmp_z, df_fin$prcp_z)
cor(df_fin$tmp, df_fin$prcp)

plot(df_fin$tmp_z, df_fin$prcp_z)
plot(df_fin$tmp, df_fin$prcp)



summary(gam_model)
plot(gam_model, page = 1, shade = T)
appraise(gam_model)
concurvity(gam_model)
gam.check(gam_model)
k.check(gam_model)

concurvity(gam_model)

# FRAMEWORKD: stepwise preicor selection -------------------------------------------------------------

# Load necessary libraries
library(mgcv)
library(ggplot2)

# filter the extreme values: 
# Filter out extreme values or cap them
df_fin_filtered <- df_fin[df_fin$sum_stems_juvenile < quantile(df_fin$sum_stems_juvenile, 0.99), ]




# Step 1: Define the base model with spatial smoothers (ti(x, y))
base_model <- gam(sum_stems_juvenile ~ ti(x, y), 
                  family = Tweedie(p = 1.46), 
                  method = 'REML',
                  data = df_fin_filtered)
summary(base_model)

# Step 2: Add temperature smoother
temp_model <- gam(sum_stems_juvenile ~ ti(x, y) + s(tmp_z, k = 7), 
                  family = Tweedie(p = 1.46), 
                  method = 'REML',
                  data = df_fin_filtered)
summary(temp_model)

# Step 3: Add precipitation smoother
precip_model <- gam(sum_stems_juvenile ~ ti(x, y) + s(prcp_z, k = 7), 
                    family = Tweedie(p = 1.46), 
                    method = 'REML',
                    data = df_fin_filtered)
summary(precip_model)

# Step 4: Combine temperature and precipitation smoothers
temp_precip_model <- gam(sum_stems_juvenile ~ ti(x, y) + s(tmp_z, k = 7) + s(prcp_z, k = 7), 
                         family = Tweedie(p = 1.46), 
                         method = 'REML',
                         data = df_fin_filtered)
summary(temp_precip_model)

# Step 5: Add interaction term between temperature and precipitation
interaction_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10), 
                         family = Tweedie(p = 1.46), 
                         method = 'REML',
                         data = df_fin_filtered)
summary(interaction_model)

# Step 6: Add random effects for country and region
random_effects_model <- gam(sum_stems_juvenile ~ ti(x, y) + s(tmp_z, k = 7) + s(prcp_z, k = 7) +
                              s(country_abbr, bs = "re") + s(region, bs = "re"),
                            family = Tweedie(p = 1.46), 
                            method = 'REML',
                            data = df_fin_filtered)
summary(random_effects_model)

# Step 7: Add management intensity by country
management_model <- gam(sum_stems_juvenile ~ ti(x, y) + s(tmp_z, k = 7) + s(prcp_z, k = 7) +
                          s(country_abbr, bs = "re") + s(region, bs = "re") + 
                          s(management_intensity, by = country_pooled, k = 10),
                        family = Tweedie(p = 1.46), 
                        method = 'REML',
                        data = df_fin_filtered)
summary(management_model)

# Step 8: Compare models using AIC
AIC(base_model, temp_model, precip_model, temp_precip_model, interaction_model, random_effects_model, management_model)

# Step 9: Check residual plots and diagnostics for the final model
best_model <- interaction_model      # Choose your best model after comparison
appraise(best_model)
plot(best_model, page=1)

# Optional: You can check concurvity or k.check if needed
concurvity(best_model)
k.check(best_model)


# continue with teh improved models: ------------------

# Load necessary libraries
library(mgcv)
library(ggplot2)

# Step 1: Base model with interaction between spatial variables (x, y) and temp/precip interaction
base_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10), 
                  family = Tweedie(p = 1.46), 
                  data = df_fin_filtered)  # Use filtered data to remove extreme values
summary(base_model)

# Step 2: Add tmp_z as an individual smoother
temp_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10) + 
                    s(tmp_z, k = 7), 
                  family = Tweedie(p = 1.46), 
                  data = df_fin_filtered)
summary(temp_model)

# Step 3: Add prcp_z as an individual smoother
precip_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10) + 
                      s(prcp_z, k = 7), 
                    family = Tweedie(p = 1.46), 
                    data = df_fin_filtered)
summary(precip_model)

# Step 4: Add both tmp_z and prcp_z as individual smoothers
temp_precip_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10) + 
                           s(tmp_z, k = 7) + s(prcp_z, k = 7), 
                         family = Tweedie(p = 1.46), 
                         data = df_fin_filtered)
summary(temp_precip_model)

# Step 5: Add random effects for region and country
random_effects_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10) + 
                              s(tmp_z, k = 7) + s(prcp_z, k = 7) + 
                              s(region, bs = 're') + s(country_abbr, bs = 're'), 
                            family = Tweedie(p = 1.46), 
                            data = df_fin_filtered)
summary(random_effects_model)

# Step 6: Add management intensity as a predictor (considering interaction with country)
management_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10) + 
                          s(tmp_z, k = 7) + s(prcp_z, k = 7) + 
                          s(region, bs = 're') + s(country_abbr, bs = 're') + 
                          s(management_intensity, by = country_pooled, k = 10), 
                        family = Tweedie(p = 1.46), 
                        data = df_fin_filtered)
summary(management_model)

# Step 7: Add additional predictors step by step
# Adding distance_edge
distance_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10) + 
                        s(tmp_z, k = 7) + s(prcp_z, k = 7) + 
                        s(region, bs = 're') + s(country_abbr, bs = 're') + 
                        s(management_intensity, by = country_pooled, k = 10) + 
                        s(distance_edge, k = 10), 
                      family = Tweedie(p = 1.46), 
                      data = df_fin_filtered)
summary(distance_model)

# Adding disturbance_severity
disturbance_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10) + 
                           s(tmp_z, k = 7) + s(prcp_z, k = 7) + 
                           s(region, bs = 're') + s(country_abbr, bs = 're') + 
                           s(management_intensity, by = country_pooled, k = 10) + 
                           s(distance_edge, k = 10) + 
                           s(disturbance_severity, k = 10), 
                         family = Tweedie(p = 1.46), 
                         data = df_fin_filtered)
summary(disturbance_model)

# Adding sand_extract
sand_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10) + 
                    s(tmp_z, k = 7) + s(prcp_z, k = 7) + 
                    s(region, bs = 're') + s(country_abbr, bs = 're') + 
                    s(management_intensity, by = country_pooled, k = 10) + 
                    s(distance_edge, k = 10) + 
                    s(disturbance_severity, k = 10) + 
                    s(sand_extract, k = 10), 
                  family = Tweedie(p = 1.46), 
                  data = df_fin_filtered)
summary(sand_model)

# Adding clay_extract
clay_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10) + 
                    s(tmp_z, k = 7) + s(prcp_z, k = 7) + 
                    s(region, bs = 're') + s(country_abbr, bs = 're') + 
                    s(management_intensity, by = country_pooled, k = 10) + 
                    s(distance_edge, k = 10) + 
                    s(disturbance_severity, k = 10) + 
                    s(sand_extract, k = 10) + 
                    s(clay_extract, k = 10), 
                  family = Tweedie(p = 1.46), 
                  data = df_fin_filtered)
summary(clay_model)

# Adding av.nitro
nitro_model <- gam(sum_stems_juvenile ~ ti(x, y) + ti(tmp_z, prcp_z, k = 10) + 
                     s(tmp_z, k = 7) + s(prcp_z, k = 7) + 
                     s(region, bs = 're') + s(country_abbr, bs = 're') + 
                     s(management_intensity, by = country_pooled, k = 10) + 
                     s(distance_edge, k = 10) + 
                     s(disturbance_severity, k = 10) + 
                     s(sand_extract, k = 10) + 
                     s(clay_extract, k = 10) + 
                     s(av.nitro, k = 10), 
                   family = Tweedie(p = 1.46), 
                   data = df_fin_filtered)
summary(nitro_model)

# Step 8: Compare AIC values after adding each predictor
AIC(base_model, temp_model, precip_model, temp_precip_model, interaction_model, 
    random_effects_model, management_model, distance_model, disturbance_model, 
    sand_model, clay_model, nitro_model)

# Step 9: Perform diagnostics on the best model based on AIC
best_model <- nitro_model  # Replace with the best model after AIC comparison
appraise(best_model)  # Model diagnostics
plot(best_model)      # Plot smoothers

# Optional: Check concurvity and k-check for potential collinearity and overfitting
concurvity(best_model)
k.check(best_model)






# END

appraise(model_nb)
summary(model_nb)
plot(model_nb, page = 1, shade = T)

# Use dredge to rank all models based on AIC
model_set <- dredge(global_model_juv)

# View the top models
head(model_set)

# Select the best model based on AIC
best_model <- get.models(model_set, 2)[[1]]

# Summary of the best model
summary(best_model)


appraise(best_model)
summary(best_model)
plot(best_model, page = 1, shade = T)
gam.check(best_model)
k.check(best_model)



# quick visualization -----------------------------------------------------------
##### quick plotting -----------------------------------------------
fin.m <- best_model

# check for tempoeal autocorrelation

# Extract residuals from the model
#residuals <- resid(m.agg.previous$lme, type = "normalized")

# Plot ACF of the residuals
#acf(residuals, main="ACF of Model Residuals")
#pacf(residuals, main="ACF of Model Residuals")


# plot in easy way how tdoes the k value affect interaction
p1 <- ggpredict(fin.m, terms = "drought_spei12 [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m, terms = "tmp_z [all]", allow.new.levels = TRUE)
p3 <- ggpredict(fin.m, terms = "prcp_z [all]", allow.new.levels = TRUE)
#p4 <- ggpredict(fin.m, terms = c("tmp_z", "prcp_z [-1, 0, 1]"), allow.new.levels = TRUE)
p4 <- ggpredict(fin.m, terms = c("tmp_z", "prcp_z"), allow.new.levels = TRUE)


#p_df <- as.data.frame(p3)
# test simple plot:
plot1<-ggplot(p1, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  theme_classic2()

plot2<-ggplot(p2, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  theme_classic2()

plot3<-ggplot(p3, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  theme_classic2()

plot4<-ggplot(p4, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  geom_line(aes(color = group, linetype = group), linewidth = 1) +
  theme_classic2()

(plot4)
ggarrange(plot1,
          plot2,plot3,
          #plot4, 
          ncol = 2, nrow = 2)


# 3D interaction plot for temperature and precipitation
vis.gam(best_model, view = c("tmp_z", "prcp_z"), plot.type = "persp", theta = 30, phi = 30)




### drivers saplings -----------------------------------------------------------------------------
# Fit a global model with all possible predictors
global_model_sapl <- gam(sum_stems_sapling ~ #s(drought_spei6, k = 8) +
                           s(tmp_z, k = 7) +
                           s(prcp_z, k = 7) +
                           # s(drought_spei12, k = 10) +
                           #s(drought_spei24, k = 10) +
                           #  s(distance_edge, k = 10) +  # did not improved the mode
                           disturbance_severity +
                           sand_extract  +
                           s(clay_extract, k = 10) +
                           s(depth_extract, k = 10) +
                           s(av.nitro, k = 15) +
                           s(country_abbr, bs = "re"),
                         family = tw(), 
                         data = df_fin,
                         method = "REML",
                         na.action = na.fail)


appraise(global_model_sapl)
summary(global_model_sapl)
plot(global_model_sapl, page = 1, shade = T)

# Use dredge to rank all models based on AIC
model_set <- dredge(global_model_sapl)

# View the top models
head(model_set)

# Select the best model based on AIC
best_model <- get.models(model_set, 2)[[1]]

# Summary of the best model
summary(best_model)


appraise(best_model)
summary(best_model)
plot(best_model, page = 1, shade = T)
gam.check(best_model)
k.check(best_model)














## OLD ------------

df_model <- df_fin %>% 
  mutate(stem_density = as.integer(stem_density)) %>% 
  # create a factor for country based on region ID
  mutate(
    region = str_extract(cluster, "^\\d{2}"), # Extracts the first two digits
    country = case_when(
      region %in% c(11, 12, 14, 18, 19, 20, 25) ~ 'DE', # Germany
      region == 17 ~ 'PL', # Poland
      region %in% c(15, 26) ~ 'CZ', # Czech Republic
      region == 13 ~ 'AT', # Austria
      region == 16 ~ 'SK', # Slovakia
      region == 23 ~ 'SI', # Slovenia
      region == 21 ~ 'IT', # Italy
      region == 22 ~ 'CH', # Switzerland
      region %in% c(24, 27) ~ 'FR', # France
      TRUE ~ NA_character_ # For any region not specified above
    )
  ) %>% 
  mutate(country = factor(country),
         dominant_species = factor(dominant_species)
         #region = as.character(region)
  )

# see differences among countries
df_model %>% 
  ggplot(aes(y = stem_density,
             x = country)) + 
  geom_boxplot()

str(df_model)

# test predictors
ggplot(df_model, aes(x = tmp, y = stem_density)) +
  geom_point() +
  geom_smooth(method = "loess", se = F) +
  theme_bw()


# test glm with neg bin family: put all predictors together
library(MASS)

# which one is better: simple management intensity or salvage intensity?
m0.manag <-glm.nb(stem_density ~ management_intensity, data = df_model, na.action = "na.fail") # + tmp_z 
m0.salvage <-glm.nb(stem_density ~ salvage_intensity , data = df_model, na.action = "na.fail") # + tmp_z 
m0.protect <-glm.nb(stem_density ~ protection_intensity, data = df_model, na.action = "na.fail") # + tmp_z 

# add country effect
m1.manag <-glm.nb(stem_density ~ management_intensity + country, data = df_model, na.action = "na.fail") # + tmp_z 
m1.salvage <-glm.nb(stem_density ~ salvage_intensity + country, data = df_model, na.action = "na.fail") # + tmp_z 
m1.protect <-glm.nb(stem_density ~ protection_intensity + country, data = df_model, na.action = "na.fail") # + tmp_z 

m2 <- glm.nb(stem_density ~ salvage_intensity + protection_intensity, data = df_model, na.action = "na.fail")
m2.country <- glm.nb(stem_density ~ salvage_intensity + protection_intensity + country, data = df_model, na.action = "na.fail")

m.country <- glm.nb(stem_density ~ country, data = df_model, na.action = "na.fail")

AICc(m0.manag, m0.salvage, m0.protect, m1.manag, m1.salvage, m1.protect,m.country, m2.country, m2)
BIC(m0.manag, m0.salvage, m0.protect, m1.manag, m1.salvage, m1.protect)

# seemms that simple managemnet intensity is better then salvage and protection intensity


plot(allEffects(m1.manag))

# put together all of predictors to decide which one is teh best: use only management_intensity
# kep only climate predictors that are less correlated (cor coeff <0.3) -> remove SPEI

global.negbin1 <- glm.nb(stem_density ~ tmp  + tmp_z + prcp_z + prec +
                           disturbance_severity + distance_edge  + 
                           management_intensity + 
                           clay_extract + 
                           country+
                           sand_extract + depth_extract + av.nitro, data = df_model, na.action = "na.fail") # + tmp_z 

# remove the largest VIF: kep only z-scores! 
global.negbin2 <- glm.nb(stem_density ~ tmp_z + 
                           prcp_z+ #spei + prec +
                           disturbance_severity + distance_edge  + 
                           management_intensity + 
                           clay_extract + 
                           #salvage_intensity + protection_intensity + 
                           country+
                           sand_extract + depth_extract + av.nitro, data = df_model, na.action = "na.fail") # + tmp_z 

vif(global.negbin2)
# remove non significant precdictors: sannd, depth, av.nitro
global.negbin3 <- glm.nb(stem_density ~ tmp_z + 
                           prcp_z+ #spei + prec +
                           disturbance_severity + distance_edge  + 
                           management_intensity + 
                           clay_extract + 
                           #salvage_intensity + protection_intensity + 
                           country #+
                         #  sand_extract + depth_extract + av.nitro
                         , data = df_model, na.action = "na.fail")

vif(global.negbin3)

# add inuteraction ebtween temp and spei, remove prec
global.negbin4 <- glm.nb(stem_density ~ tmp *spei +
                           prcp_z +# spei + prec +  
                           disturbance_severity + distance_edge  + 
                           management_intensity + 
                           clay_extract + 
                           #salvage_intensity + protection_intensity + 
                           country, #+
                         #  sand_extract + depth_extract + av.nitro
                         data = df_model, na.action = "na.fail") # + tmp_z 

# remove spei, highly correlayted tiwth tmp and prec

global.negbin5 <- glm.nb(stem_density ~ tmp  + #tmp_z + 
                           prcp_z+ 
                           #spei + 
                           prec +
                           disturbance_severity + distance_edge  + 
                           management_intensity + 
                           clay_extract + 
                           #salvage_intensity + protection_intensity + 
                           country, #+
                         #  sand_extract + depth_extract + av.nitro
                         data = df_model, na.action = "na.fail") # + tmp_z 

# keep only z-scores
global.negbin6 <- glm.nb(stem_density ~ tmp_z  + #tmp_z + 
                           prcp_z+ 
                           #spei + 
                           #prec +
                           disturbance_severity + distance_edge  + 
                           management_intensity + 
                           clay_extract + 
                           #salvage_intensity + protection_intensity + 
                           country, #+
                         #  sand_extract + depth_extract + av.nitro
                         data = df_model, na.action = "na.fail") # + tmp_z 







vif(global.negbin6)

AIC(global.negbin1,global.negbin2,global.negbin3,global.negbin4, global.negbin5, global.negbin6)

summary(global.negbin5)
vif(global.negbin2)
simulateResiduals(global.negbin3, plot = T)

#vif(global.negbin)
glob1 <- dredge(global.negbin)



# best predictors are: clay, distance_edge. dist_severity, tmp_z keep also managemebt and countries. add 

# tmp_ is a better predictors than spei alone
global.negbin2 <- glm.nb(stem_density ~  tmp_z + #prcp_z + 
                           disturbance_severity + distance_edge  + 
                           management_intensity + clay_extract , data = df_model, na.action = "na.fail") # + tmp_z 


dredge(global.negbin2)

vif(global.negbin2)

plot(allEffects(global.negbin2))
summary(global.negbin2)

simulateResiduals(global.negbin2, plot = T)




# include random effects of countries:
library(glmmTMB)

glmm1 <- glmmTMB(stem_density ~ tmp  + tmp_z + prcp_z + prec + spei +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   country+
                   sand_extract + depth_extract + av.nitro +
                   (1 | country), data = df_model, 
                 family = nbinom2, na.action = "na.fail")

glmm2 <- glmmTMB(stem_density ~ tmp*spei + prcp_z +# prec +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   # country+
                   # sand_extract + depth_extract + av.nitro +
                   (1 | country), data = df_model, 
                 family = nbinom2, na.action = "na.fail")

glmm3 <- glmmTMB(stem_density ~ tmp*spei + prcp_z +# prec +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   # country+
                   # sand_extract + depth_extract + av.nitro +
                   (1 | country/management_intensity), data = df_model, 
                 family = nbinom2, na.action = "na.fail")

# use tmp_z as well
glmm4 <- glmmTMB(stem_density ~ tmp_z*spei + prcp_z +# prec +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   # country+
                   # sand_extract + depth_extract + av.nitro +
                   (1 | country), data = df_model, 
                 family = nbinom2, na.action = "na.fail")


testZeroInflation(glmm4)
# the p is significant, so it has too many zeros!


# use zero inflated model: ----------------------------------------------------

# check which one of clim data are the best predistors:
m.zi.clim <- glmmTMB(stem_density ~ tmp  + tmp_z + prcp_z + prec + spei + #tmp:spei +
                       # disturbance_severity + distance_edge  + 
                       # # management_intensity + 
                       # clay_extract + 
                       # #country+
                       # sand_extract + depth_extract + av.nitro +
                       (1 | country/management_intensity), 
                     ziformula = ~1,
                     data = df_model, 
                     family = nbinom2, na.action = "na.fail")


zi.clim.full <- dredge(m.zi.clim)

m.zi0 <- glmmTMB(stem_density ~ tmp  + #tmp_z + prcp_z + 
                   prec + spei + #tmp:spei +
                   disturbance_severity + distance_edge  + 
                   # management_intensity + 
                   clay_extract + 
                   #country+
                   sand_extract + depth_extract + av.nitro +
                   (1 | country/management_intensity), 
                 ziformula = ~1,
                 data = df_model, 
                 family = nbinom2, na.action = "na.fail")

summary(m.zi0)

m.zi01 <- glmmTMB(stem_density ~ tmp  + #tmp_z + prcp_z + 
                    # prec + 
                    spei + #tmp:spei +
                    disturbance_severity + distance_edge  + 
                    # management_intensity + 
                    clay_extract + 
                    #country+
                    #sand_extract + depth_extract + av.nitro +
                    (1 | country/management_intensity), 
                  ziformula = ~1,
                  data = df_model, 
                  family = nbinom2, na.action = "na.fail")

AIC(m.zi0, m.zi01)

m.zi01 <- glmmTMB(stem_density ~ tmp  + #tmp_z + prcp_z + 
                    prec + 
                    spei + #tmp:spei +
                    disturbance_severity + distance_edge  + 
                    # management_intensity + 
                    clay_extract + 
                    #country+
                    #sand_extract + depth_extract + av.nitro +
                    (1 | country/management_intensity), 
                  ziformula = ~1,
                  data = df_model, 
                  family = nbinom2, na.action = "na.fail")

fin.m <- m.zi01  # the best, also has good residuals
# add poly, only meaningful for prec
m.zi02 <- glmmTMB(stem_density ~ poly(tmp,2)  + #tmp_z + prcp_z + 
                    poly(prec,2) + 
                    poly(spei,2) + #tmp:spei +
                    disturbance_severity + distance_edge  + 
                    # management_intensity + 
                    clay_extract + 
                    #country+
                    #sand_extract + depth_extract + av.nitro +
                    (1 | country/management_intensity), 
                  ziformula = ~1,
                  data = df_model, 
                  family = nbinom2, na.action = "na.fail")

# add poly, only meaningful for prec - cnvergece issues, remove
m.zi03 <- glmmTMB(stem_density ~ tmp  + #tmp_z + prcp_z + 
                    poly(prec,2) + 
                    spei + #tmp:spei +
                    disturbance_severity + distance_edge  + 
                    # management_intensity + 
                    clay_extract + 
                    #country+
                    #sand_extract + depth_extract + av.nitro +
                    (1 | country/management_intensity), 
                  ziformula = ~1,
                  data = df_model, 
                  family = nbinom2, na.action = "na.fail")


AIC(m.zi02, m.zi01,m.zi03)
summary(m.zi01)
plot(allEffects(m.zi01))

simulateResiduals(m.zi01, plot = T)

library(sjPlot)
sjPlot::tab_model(m.zi01,
                  #col.header = c(as.character(qntils), 'mean'),
                  #show.rownames = F,
                  file="outTable/drivers.doc") 


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

