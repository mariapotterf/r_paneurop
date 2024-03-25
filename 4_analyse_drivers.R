# create final database
# merge veg data
# predictors - climate, disturbance chars, soils
# dependent - regeneretion, density, richness
# output file
# 

# Libs --------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(terra)
library(ggplot2)
library(lubridate)
library(stringr)

library(ggpubr)

library(dunn.test) # Duncan test, post hoc Kruskal-Wallis


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

library(AER)  # for hurdle model


# read data --------------------------------------------------------------------

climate           <- fread("outData/climate_18_23.csv")
spei              <- fread("outData/SPEI.csv")
distance_to_edge  <- fread("outData/distance_to_edge.csv")
disturbance_chars <- fread("outData/disturbance_chars.csv")
soil              <- as.data.frame(terra::vect("rawData/extracted_soil_data/extracted_soil_data_completed.gpkg"))
terrain           <- fread("outData/terrain.csv")

# get vegetation data
load("outData/veg.Rdata")


# clean up data -----------------------------------------------------------------
# sum stems first up across vertical groups!!

# get only samplings
df_stems_sapl <- stem_dens_ha_cluster_sum %>% 
  filter(VegType == 'Regeneration') %>% 
  group_by(cluster, management_intensity) %>% 
  summarise(sum_stems = sum(total_stems))# %>% 

# get saplings and juveniles
df_stems_both <- stem_dens_ha_cluster_sum %>% 
  filter(VegType != 'Survivor') %>% 
  group_by(cluster, management_intensity) %>% 
  summarise(sum_stems = sum(total_stems))# %>% 

table(stem_dens_species_long$VegType,stem_dens_species_long$Species)


# aggregate by cluster: sum up teh species, by height category
stem_dens_species_cluster_long_share <- stem_dens_species_long %>% 
  ungroup(.) %>% 
 # dplyr::select(-management_intensity) %>% 
  group_by(cluster, Species, VegType, management_intensity) %>%
    summarise(stem_density = sum(stem_density))
    

#fwrite(stem_dens_species_cluster_long_share, 'vegData_cluster_densities.csv')

disturbance_chars <- disturbance_chars %>%
  mutate(country = as.character(country))

#left_join(dplyr::filter(spei, year %in% 2018:2023)) #%>% 

#Process SPEI:  get cluster number
spei <- spei %>% 
  dplyr::mutate(year  = lubridate::year(date),
                month = lubridate::month(date)) %>%  
  dplyr::select(-c(date, -scale)) %>% 
  mutate(cluster = str_sub(ID, 4, -3))#%>%


spei_ID <- spei %>% 
  ungroup(.) %>% 
  dplyr::filter(year %in% 2018:2023 ) %>% # & month %in% 4:9# select just vegetation season
  dplyr::filter_all(all_vars(!is.infinite(.))) %>% # remove all infinite values
  #View()
  group_by(ID, year) %>%
  summarise(spei = mean(spei, na.rm = T)) #%>% 


spei_cluster <- spei %>% 
  ungroup(.) %>% 
  dplyr::filter(year %in% 2018:2023 ) %>% # & month %in% 4:9# select just vegetation season
  dplyr::filter_all(all_vars(!is.infinite(.))) %>% # remove all infinite values
  #View()
  group_by(cluster, year) %>%
  summarise(spei = mean(spei, na.rm = T)) #%>% 



# merge all predictors together, ID level
df_predictors <- 
  climate %>% 
  left_join(dplyr::select(distance_to_edge, c(-country, -dist_ID ))) %>% 
  left_join(dplyr::select(disturbance_chars, c(-region, -country ))) %>% 
  left_join(soil) %>%
  left_join(spei_ID) %>% 
  left_join(dplyr::select(terrain, c(-country, -region, -cluster.x, -cluster.y))) #%>% 
  #mutate(cluster = str_sub(ID, 4, -3)) 

length(unique(df_predictors$cluster)) # 957
anyNA(df_predictors)
# merge predcitors on cluster level: calculate means
# keep only temp and prcp: 2021-2023 average
df_predictors_cluster <- 
  df_predictors %>% 
  dplyr::filter(year %in% 2018:2023) %>% 
  group_by(cluster) %>% 
  summarise(tmp = mean(tmp, na.rm = T),
            prec = mean(prec, na.rm = T),
            tmp_z = mean(tmp_z, na.rm = T),
            prcp_z = mean(prcp_z, na.rm = T),
            spei   = mean(spei, na.rm = T),
            distance = mean(distance, na.rm = T),
            disturbance_year= as.integer(mean(disturbance_year, na.rm = T)),
            disturbance_severity= mean(disturbance_severity, na.rm = T),
            disturbance_agent= as.integer(mean(disturbance_agent, na.rm = T)),
            elevation = mean(elevation, na.rm = T),
            #region = mean(tmp, na.rm = T),
            #country = mean(tmp, na.rm = T),
            sand_extract= mean(sand_extract, na.rm = T),
            silt_extract = mean(silt_extract, na.rm = T),
            clay_extract = mean(clay_extract, na.rm = T),
            depth_extract = mean(depth_extract, na.rm = T),
            av.nitro    = mean(av.nitro, na.rm = T),
            slope   = mean(slope, na.rm = T),
            aspect= mean(aspect, na.rm = T))


df_predictors_ID <- 
  df_predictors %>% 
  dplyr::filter(year %in% 2018:2023) %>% 
  group_by(ID) %>% 
  summarise(tmp = mean(tmp, na.rm = T),
            prec = mean(prec, na.rm = T),
            tmp_z = mean(tmp_z, na.rm = T),
            prcp_z = mean(prcp_z, na.rm = T),
            spei   = mean(spei, na.rm = T),
            distance = mean(distance, na.rm = T),
            disturbance_year= as.integer(mean(disturbance_year, na.rm = T)),
            disturbance_severity= mean(disturbance_severity, na.rm = T),
            disturbance_agent= as.integer(mean(disturbance_agent, na.rm = T)),
            elevation = mean(elevation, na.rm = T),
            #region = mean(tmp, na.rm = T),
            #country = mean(tmp, na.rm = T),
            sand_extract= mean(sand_extract, na.rm = T),
            silt_extract = mean(silt_extract, na.rm = T),
            clay_extract = mean(clay_extract, na.rm = T),
            depth_extract = mean(depth_extract, na.rm = T),
            av.nitro    = mean(av.nitro, na.rm = T),
            slope   = mean(slope, na.rm = T),
            aspect= mean(aspect, na.rm = T))


fwrite(df_predictors_ID, 'outData/all_predictors_ID.csv')

#fwrite(df_predictors_cluster, 'outData/all_predictors_cluster.csv')

# merge with vegetation data
# select the dominant species per cluster
df_IVI_max <- df_IVI %>% 
  dplyr::select(cluster, Species, country,rIVI) %>% 
  ungroup(.) %>% 
  group_by(cluster, country) %>% 
  filter(rIVI == max(rIVI)) %>% 
  slice(1)

# cluster level
df_cluster_veg <-
  df_IVI_max %>% 
  left_join(df_richness) %>%
  left_join(df_vert) %>%
  left_join(df_stems)# %>%

anyNA(df_cluster_veg)

# cluster data share with predictors:
# rename not necessary predictors for clustering, rename veg variables
df_cluster_full <- df_cluster_veg %>% 
  left_join(df_predictors_cluster, by = join_by(cluster)) %>% 
  ungroup() %>% 
  dplyr::select(-c(disturbance_year, 
                   disturbance_agent,
                   silt_extract, country,
                   elevation,
                   slope, aspect)) %>% 
  dplyr::rename(dominant_species = Species, 
                stem_density = sum_stems,
                distance_edge = distance,
                n_vertical = n_layers)
  

fwrite(df_cluster_full, 'outData/indicators_for_cluster_analysis.csv')

length(unique(df_cluster$cluster))
length(unique(df_cluster_full$cluster)) # 849!   - final clusters, 4-5 plots
length(unique(df_predictors_cluster$cluster))  # 957 - all clusters, from even with less plots


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
# higher management intensity, higher species richness
# higher managemnet, less of teh advanced regeneration - 
#   due to the salvage logging and following 

# drivers:
# apply drivers for all 4 veg indicators??
# effect of climate: 
#   - higher drought, less regeneration; higher change of delayed regeneration

# effect of disturbance size:
#   - higher patch size, less regeneration
#   - higher severity, low regenerations

# environmental condistins;
#  soil: 
#    - more stem density at better soil conditions (higher depth, nutrient content)


# Make 2d denisty plots: 

# Merge structure & composition into single space 


View(df_cluster_full)
# get summary for the scatter plot, one point is one country & management
df_summary <- df_cluster_full %>%
 # group_by(manag) %>% # country, 
  summarize(
    rich_mean = median(rIVI, na.rm = TRUE),
    rich_sd = sd(rIVI, na.rm = TRUE),
    rich_min = quantile(rIVI, 0.25, na.rm = TRUE), # rich_mean -rich_sd, #
    rich_max = quantile(rIVI, 0.75, na.rm = TRUE), # rich_mean +rich_sd,
    dens_mean = median(stem_density   , na.rm = TRUE),
    dens_sd   = sd(stem_density  , na.rm = TRUE),
    dens_min = quantile(stem_density  , 0.25, na.rm = TRUE), # dens_mean - dens_sd, #
    dens_max = quantile(stem_density  , 0.75, na.rm = TRUE), # dens_mean + dens_sd, #
    .groups = 'drop'
  )

(df_summary)


# test general trends of vegetation along management gradient: -------------------------------
p1 <- df_cluster_full %>% 
  ggplot(aes(x = management_intensity ,
             y = rIVI))+ 
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = 'loess')  
  
  
p2 <- df_cluster_full %>% 
    ggplot(aes(x = management_intensity,
               y = richness))+ 
  geom_jitter(alpha = 0.5) +
    geom_smooth()

p3 <- df_cluster_full %>% 
  ggplot(aes(x = management_intensity,
             y = stem_density))+ 
  geom_jitter(alpha = 0.5) +
  geom_smooth()

p4 <- df_cluster_full %>% 
  ggplot(aes(x = management_intensity,
             y = n_vertical))+
  geom_jitter(alpha = 0.5) +
  geom_smooth()
  


ggarrange(p1, p2, p3,p4, nrow = 2, ncol = 2)


# try simple glm: values between 0-1 - use beta distribution:
# # make sure that no dat are equal 0 or 1
# df_cluster_full
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




# prepare plot with raster -------------------------------

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



# Define the function to prepare density data --------------------------
prepare_density_data <- function(df_cluster_full, 
                                 x = 'rIVI', y = 'stem_density',
                                 percentages = c(0.25, 0.50, 0.75, 0.90, 1)) {
  # Calculate 2D density
  d <- kde2d(df_cluster_full[[x]], df_cluster_full[[y]], n = 300)
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

# Density plots: each variable against another one --------------------
# 1. rIVI vs stem_density
# 2. rIVI vs richness
# 3. rIVI vs n_vertical
# 4. stem_density vs richness
# 5. stem_density vs n_vertical
# 6. richness vs n_vertical 


### rIVI vs sum_stems -----------------------------------------
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     x = 'stem_density',
                                     y = 'rIVI')


p_IVI_stems <- make_2D_plot(df_plot_data, 
             x_lab = "Stem density [n/ha]",
             y_lab = "Importance value [%]")

(p_IVI_stems)
### rIVI vs richnes -----------------------------------
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     x = 'rIVI', 
                                     y = 'richness')


p_IVI_richness <- make_2D_plot(df_plot_data, 
                            x_lab = "Importance value [%]", # Default axis labels, can be customized
                            y_lab = "Richness [counts]"  )


### rIVI vs vertical -----------------------------------
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     x = 'rIVI', 
                                     y = 'n_vertical')


p_IVI_vert <- make_2D_plot(df_plot_data, 
                               x_lab = "Importance value [%]", # Default axis labels, can be customized
                               y_lab = "# vertical layers [counts]"  )



### stem_density vs richness -----------------------------------
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     x = 'stem_density', 
                                     y = 'richness')


p_dens_richness <- make_2D_plot(df_plot_data, 
                           x_lab = "Stem density [n/ha]", # Default axis labels, can be customized
                           y_lab = "Richness [counts]"  )



### stem_density vs n_vertcal -----------------------------------
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     x = 'stem_density', 
                                     y = 'n_vertical')


p_dens_vertical <- make_2D_plot(df_plot_data, 
                                x_lab = "Stem density [n/ha]", # Default axis labels, can be customized
                                y_lab = "# vertical layers [counts]"  )


### n_vertical vs richness -----------------------------------
df_plot_data <- prepare_density_data(df_cluster_full, 
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
          p_vertical_richness, nrow = 2, ncol = 3, common.legend = T, legend = 'bottom', align = 'hv')

# DEnsity plot wth management ---------------------------------------------
###### rIVI vs management ------------------------------------------------
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     y = 'rIVI', 
                                     x = 'management_intensity')


p_IVI_manag <- make_2D_plot(df_plot_data, 
                          # x = 'x', # Default column names, can be changed
                          # y = 'y', # Default column names, can be changed
                           y_lab = "Importance value [%]", # Default axis labels, can be customized
                           x_lab = "Management intensity [%]"  )

(p_IVI_manag)

###### richness vs management ------------------------------------------------
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     y = 'richness', 
                                     x = 'management_intensity')


p_richness_manag <- make_2D_plot(df_plot_data, 
                           # x = 'x', # Default column names, can be changed
                           # y = 'y', # Default column names, can be changed
                            y_lab = "Richness [#]", # Default axis labels, can be customized
                            x_lab = "Management intensity [%]"  )

(p_richness_manag)



###### stem_density vs management ------------------------------------------------
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     y = 'stem_density', 
                                     x = 'management_intensity')


p_stem_dens_manag <- make_2D_plot(df_plot_data, 
                                # x = 'x', # Default column names, can be changed
                                # y = 'y', # Default column names, can be changed
                                 y_lab = "Stem density [n/ha]", # Default axis labels, can be customized
                                 x_lab = "Management intensity [%]"  )

(p_stem_dens_manag)




###### n_vertical vs management ------------------------------------------------
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     y = 'n_vertical', 
                                     x = 'management_intensity')


p_vertical_manag <- make_2D_plot(df_plot_data, 
                                 # x = 'x', # Default column names, can be changed
                                 # y = 'y', # Default column names, can be changed
                                  y_lab = "Vertical layers [#]", # Default axis labels, can be customized
                                  x_lab = "Management intensity [%]"  )

(p_vertical_manag)



# Merge plots indicators vs management intensity --------------------------

ggarrange(p_IVI_manag, p_richness_manag, 
          p_stem_dens_manag, p_vertical_manag, 
          nrow = 2, ncol = 2,
          common.legend = T, #legend.position = 'right', 
          align = 'hv')




# get summary table per quantiles -------------------------------------------------

# Summary table  ----------------------------------------------------------

# Represent results using quantiles, as they are skewed?
qntils = c(0, 0.25, 0.5, 0.75, 1)

qs_dat_fin <- 
  df_cluster_full %>% 
  ungroup(.) %>% 
  #filter(year %in% 2015:2021) %>% 
  dplyr::reframe(rIVI            = quantile(rIVI, qntils, na.rm = T ),
                 richness        = quantile(richness, qntils, na.rm = T ),
                 stem_density    = quantile(stem_density , qntils, na.rm = T ),
                 n_vertical      = quantile(n_vertical, qntils, na.rm = T )) %>% 
  t() %>%
  round(1) %>% 
  as.data.frame()

(qs_dat_fin)  


means_dat_fin <- 
  df_cluster_full %>% 
  ungroup(.) %>% 
  dplyr::reframe(rIVI            = mean(rIVI,na.rm = T ),
                 richness        = mean(richness, na.rm = T ),
                 stem_density    = mean(stem_density , na.rm = T ),
                 n_vertical      = mean(n_vertical, na.rm = T )) %>% 
  t() %>%
  round(1) %>% 
  as.data.frame() 


str(means_dat_fin)

# merge qs and mean tables
summary_out <- cbind(qs_dat_fin, means_dat_fin)# %>%
View(summary_out)

# Export as a nice table in word:
sjPlot::tab_df(summary_out,
               col.header = c(as.character(qntils), 'mean'),
               show.rownames = TRUE,
               file="outTable/summary_out.doc",
               digits = 1) 








# 2. species composition per stems: saplings vs juveniles?? --------------
# how many plots with Mature trees?
# 
dom_species <- stem_dens_species_long %>% 
  filter(VegType != 'Survivor' ) %>%  # Survivors: on 191
  group_by(Species, VegType) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup(.) %>% 
  group_by(VegType) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share = sum_stems/sum_vegType*100) #%>% 

# dominant species by country?
dom_species_country <- 
  stem_dens_species_long %>% 
  dplyr::filter(cluster %in% cluster_id_keep) %>% 
  filter(VegType != 'Survivor' ) %>%  # Survivors: on 191
  group_by(Species, VegType, manag, country ) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup(.) %>% 
  group_by(VegType,manag, country) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share = sum_stems/sum_vegType*100)



# SElect the 5 species with teh highest share per management, vertical structure
top_species_per_country <- dom_species_country %>%
  ungroup(.) %>% 
  group_by(country, VegType, manag) %>%
  arrange(desc(share)) %>%  # Arrange in descending order of share within each group
  slice_max(order_by = share, n = 5)  # Select top 5 species with max share in each group

View(top_species_per_country)



# how many species in total?
distinct_species_per_manag <- dom_species %>%
  filter(sum_stems > 0) %>%
  group_by(manag) %>%
  summarise(distinct_species_count = n_distinct(Species))

distinct_species_per_manag

# SElect the 5 species with teh highest share per management, vertical structure
top_species_per_group <- dom_species %>%
  group_by(VegType, manag) %>%
  arrange(desc(share)) %>%  # Arrange in descending order of share within each group
  slice_max(order_by = share, n = 5)  # Select top 5 species with max share in each group

(top_species_per_group)

# see all species
top_species_per_group %>% 
  #filter(
  #dplyr::filter(share > 1) %>% 
  ggplot(aes(x = manag,
             y = share,
             group = VegType,
             fill = Species)) +
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1)) +
  facet_grid(.~VegType)
  
  #geom_col('')


# Reorder Species within each group based on share
top_species_per_group <- top_species_per_group %>%
  arrange(VegType, manag, desc(share)) %>%
  mutate(Species = factor(Species, levels = unique(Species)),
         VegType = factor(VegType, 
                          levels = c("Regeneration", "advRegeneration"),
                          labels = c("Saplings", "Juveniles")))


# Create a bar plot with ordered categories
ggplot(top_species_per_group, aes(x = manag, y = share, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(share), y = share), 
            position = position_dodge(width = 0.9), 
            color = "black", size = 3, vjust = -0.5) + # Adjust vjust for label position
  facet_grid(. ~ VegType) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) +
  labs(y = "Share [%]", x = "", fill = "Species") +
  theme_classic()


ggplot(top_species_per_group, aes(x = manag, y = share, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(#aes(label = round(share,1)), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3) +
  facet_grid(. ~ VegType) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) +
  labs(y = "Share [%]", x = "Management Type", fill = "Species")


# all species by aboundance
dom_species %>% 
  dplyr::filter(share > 1) %>% 
  ggplot(aes(x = manag,
             y = share,
             fill = Species)) +
  geom_bar(position = "fill", stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#geom_col('')



# species by management
dom_species %>% 
  #dplyr::filter(VegType == 'Regeneration') %>% 
  dplyr::filter(share > 0) %>% 
  ggplot(aes(x = reorder(Species, -share),
             y = share,
             fill = VegType)) +
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(VegType~manag)
#geom_col('')


# 3. how many clusters/plots has less then 2000 regen/ha? ----------------------------------
# which country is the most affected? (where they ccur most often?)



df_stock_density <- 
  df_stems %>% 
  dplyr::filter(cluster %in% cluster_id_keep) %>% 
  mutate(dens_category = case_when(sum_stems == 0 ~ '0',
                                   sum_stems >= 1   & sum_stems < 500 ~ '1-500',
                                   sum_stems >= 500 & sum_stems < 1000 ~ '500-1000',
                                  # sum_stems >= 1000 & sum_stems < 1500 ~ '1000-1500',
                                   sum_stems >= 1000 & sum_stems < 2000 ~ '1000-2000',
                                   sum_stems >= 2000 ~ '>2000'
                                   )) %>%
  mutate(dens_category_simpl = case_when(sum_stems < 1000 ~ '<1000',
                                         sum_stems >= 1000 &sum_stems < 2000 ~ '1000-2000',
                                         sum_stems >= 2000 ~ '>2000'
  )) %>%
  mutate(dens_category = factor(dens_category, levels = 
                                  c('0', '1-500','500-1000','1000-2000', '>2000' ))) %>% 
  mutate(dens_category_simpl = factor(dens_category_simpl, levels = 
                                c('<1000','1000-2000', '>2000' ))) 
str(df_stock_density)

df_stock_density %>% 
  ggplot(aes(x = manag,
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
  filter(cluster %in% cluster_id_keep) %>%  
  dplyr::filter(stem_density>0) %>% 
  ungroup(.) %>% 
  dplyr::select(cluster, country, manag, VegType) %>%
  group_by(cluster, country, manag) %>%
  distinct(.) %>% 
  mutate(n_layers = n()) #%>% 



# Calculate frequencies of individual layers and combinations
layer_frequencies <- 
  df_vert_full %>%
  group_by(manag, cluster) %>% 
  right_join(df_stems) %>% # to account for teh empty ones
    dplyr::filter(cluster %in% cluster_id_keep) %>% 
   # View()
    mutate(VegType = case_when(
      is.na(VegType) ~ 'NO',   # Replace NA with 'NO'
      TRUE ~ VegType           # Keep existing values for non-NA cases
    )) %>% 
  summarise(layers = paste(sort(unique(VegType)), collapse = "-")) %>%
  ungroup() %>%
  count(manag, layers)


# View the result
print(layer_frequencies)

ggplot(layer_frequencies, aes(x = manag, y = n, fill = layers)) +
  geom_bar(stat = "identity", position = "fill") +
  #geom_text(aes(label = paste0(round(scales::percent, 1), "%"), 
  #              y = cumsum(n) - n/2), 
  #          position = position_fill(), 
  #          color = "white", 
  #          size = 3) +
  ylab("Percentage (%)") +
  theme_classic() +
  labs(fill = "Layers")


# 4. stems vs weather - mean 2019-2022, or SPEI ---------------------------------
# average predictors by cluster
df_predictors_clust <- df_predictors %>% 
  dplyr::filter(cluster %in% cluster_id_keep) %>% 
  group_by(cluster) %>%
  filter(year %in% 2018:2023) %>% 
  summarise(tmp                  = mean(tmp, na.rm = T),
            prec                 = mean(prec, na.rm = T), 
            tmp_z                = mean(tmp_z, na.rm = T),
            prcp_z               = mean(prcp_z, na.rm = T),
            distance             = mean(distance, na.rm = T),
            disturbance_year     = mean(disturbance_year, na.rm = T),
            disturbance_severity = mean(disturbance_severity, na.rm = T),
            disturbance_agent    = mean(disturbance_agent, na.rm = T),
            elevation            = mean(elevation, na.rm = T),
            sand_extract         = mean(sand_extract, na.rm = T),
            silt_extract         = mean(silt_extract, na.rm = T),
            clay_extract         = mean(clay_extract, na.rm = T),
            depth_extract= mean(depth_extract, na.rm = T),
            av.nitro  = mean(av.nitro, na.rm = T),
            slope = mean(slope, na.rm = T),
            aspect= mean(aspect, na.rm = T))
  
length(unique(df_predictors_clust$cluster))

  

### get average spei after disturbances ------------------------------------------
# need to update, as too many points are missing!!!
# eg ID 17_21_114_1


   
# join to stems counts: filter, as from veg data I have already excluded Italy and teh uncomplete clusters
df_fin <- #df_stems_both %>%
  df_stems %>% 
  left_join(df_predictors_clust, by = join_by(cluster)) %>% 
  left_join(spei_sub, by = join_by(cluster)) %>% 
  mutate(manag = as.factor(manag)) %>% 
  mutate(sum_stems = round(sum_stems)) %>% 
  dplyr::filter(cluster %in% cluster_id_keep) %>% 
  na.omit() %>% 
  as.data.frame()


length(unique(df_fin$cluster))

# clusters seems missing??
env_cluster <- unique(df_predictors$cluster)
veg_cluster <- unique(df_stems$cluster)

sort(env_cluster)
sort(veg_cluster)

setdiff(env_cluster, veg_cluster)

pairs(sum_stems   ~ tmp + tmp_z + distance, df_fin)

# some values are missing predictors: eg. 18_178, 18_184

#5. test model ---------------------------------

# test glm with neg bin family
library(MASS)
global.negbin <- glm.nb(sum_stems ~ tmp  + tmp_z +  spei + prec + prcp_z + disturbance_severity + distance + manag + elevation, data = df_fin, na.action = "na.fail") # + tmp_z 

dredge(global.negbin)


global.negbin2 <- glm.nb(sum_stems ~ tmp   +  prec  + spei + disturbance_severity + distance, data = df_fin, na.action = "na.fail") # + tmp_z 

dredge(global.negbin2)

vif(global.negbin2)

plot(allEffects(global.negbin2))


m.nb <- glm.nb(sum_stems ~ tmp  +  distance + manag, data = df_fin, na.action = "na.fail") # + tmp_z 

m2 <- glm.nb(sum_stems ~ tmp   +  prec  + spei + disturbance_severity + distance + manag, data = df_fin, na.action = "na.fail") # + tmp_z 

m3 <- glm.nb(sum_stems ~ tmp_z  +  distance, data = df_fin, na.action = "na.fail") # + tmp_z 

simulationOutput <- simulateResiduals(fittedModel = m2, plot = T)


cor(df_fin$tmp, df_fin$tmp_z)

AIC(m1, m2)
# zero inflated model; hurdle model:
library(pscl)
# m.zip <- zeroinfl(sum_stems ~ tmp   +  prec  + disturbance_severity + distance, #+ tmp_z 
#                   dist = "poisson", data = df_fin)

m.zinb <- zeroinfl(sum_stems ~ tmp   +  spei + prec  + disturbance_severity + distance, #+ tmp_z
                   dist = "negbin", data = df_fin)

# m.hurdle1 is the best one for now!
m.hurdle1 <- hurdle(sum_stems ~ tmp   +  spei + prec  + disturbance_severity + distance, #+ tmp_z
                   dist = "negbin", data = df_fin)

m.hurdle2 <- hurdle(sum_stems ~ tmp   +  spei + prec  + disturbance_severity + distance + manag, #+ tmp_z
                   dist = "negbin", data = df_fin)

m.hurdle3 <- hurdle(sum_stems ~  spei + disturbance_severity + distance + manag, #+ tmp_z
                    dist = "negbin", data = df_fin)

m.hurdle4 <- hurdle(sum_stems ~  spei + disturbance_severity + distance, #+ tmp_z
                    dist = "negbin", data = df_fin)

# add back temperature and precipitation, as they are from 2021-2023, spei is from 2018-2020 - the best!
m.hurdle5 <- hurdle(sum_stems ~  tmp + prec + spei + disturbance_severity + distance + manag, #+ tmp_z
                    dist = "negbin", data = df_fin)


# remove TMP as correlated with spei
m.hurdle6 <- hurdle(sum_stems ~  prec + spei + disturbance_severity + distance + manag, #+ tmp_z
                    dist = "negbin", data = df_fin)

# add soil depth
m.hurdle7 <- hurdle(sum_stems ~  tmp + prec + spei + disturbance_severity + distance + manag + depth_extract, #+ tmp_z
                    dist = "negbin", data = df_fin)

# add soil depth, remove again tmp and prec as correlated with spei
m.hurdle8 <- hurdle(sum_stems ~  spei + disturbance_severity + distance + manag + depth_extract, #+ tmp_z
                    dist = "negbin", data = df_fin)



# for several groups:
kruskal.test(sum_stems ~ manag, data = df_fin)


# Perform the Dunn test
dunnTest <- dunn.test(df_fin$sum_stems, df_fin$manag, method="bonferroni")
dunnTest



AIC(m2, m.hurdle1, m.hurdle2, m.hurdle3, m.hurdle8, m.hurdle5, m.hurdle6, m.hurdle7)

summary(m.hurdle7)
r2(m2)

# check model assumptions
plot(resid(m.hurdle1) ~ fitted(m.hurdle1))

cor(df_fin$tmp, df_fin$spei)
cor(df_fin$prec, df_fin$spei)


# go fo negative binomial, zero inf

# filtered SPEI INf values

# visualize teh hurdle model: https://library.virginia.edu/data/articles/getting-started-with-hurdle-models
sum(predict(m.hurdle1, type = "prob")[,1])  # number of zeros in teh predictoed data

# First 5 expected counts
predict(m.hurdle1, type = "response")[1:5]


#library(countreg)

#rootogram(m.hurdle1) # fit up to count xx



# make plot for hurdle5

m.hurdle5

summary(m.hurdle5)

# Function to determine if a value is significant
is_significant <- function(significance) {
  ifelse(significance == " ", "Non-significant", "Significant")
}

# Prepare data for count model
count_model_summary <- data.frame(
  Predictor = c("tmp", "prec", "spei", "disturbance_severity", "distance", "managUnmanaged"),
  Estimate = c(0.1120775, 0.0013370, 0.3777615, 0.6380436, 0.0007405, 0.1164662),
  StdError = c(0.0377988, 0.0002272, 0.3064617, 0.1849542, 0.0002089, 0.0703305),
  Significance = c("**", "***", "", "***", "***", ".")
)
count_model_summary$SigColor <- sapply(count_model_summary$Significance, is_significant)

# Plot for count model
p_count <- ggplot(count_model_summary, aes(x = Predictor, y = Estimate, color = SigColor)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - StdError, ymax = Estimate + StdError), width = 0.1) +
  scale_color_manual(values = c("Non-significant" = "grey", "Significant" = "black")) +
  geom_hline(yintercept = 0, lty = 'dashed', col = 'grey') +
  theme_classic() +
  theme(legend.position = 'none') + 
  labs(title = "Count Model Coefficients",
       y = "Estimate",
       x = "",
       color = "Significance")

# Prepare data for zero hurdle model
zero_hurdle_model_summary <- data.frame(
  Predictor = c("tmp", "prec", "spei", "disturbance_severity", "distance", "managUnmanaged"),
  Estimate = c(0.429226, 0.002572, 5.470773, -0.397865, 0.001769, -0.146980),
  StdError = c(0.194297, 0.001051, 1.794975, 0.765920, 0.001345, 0.273274),
  Significance = c("*", "*", "**", "", "", "")
)
zero_hurdle_model_summary$SigColor <- sapply(zero_hurdle_model_summary$Significance, is_significant)

# Plot for zero hurdle model
p_zero_hurdle <- ggplot(zero_hurdle_model_summary, aes(x = Predictor, y = Estimate, color = SigColor)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - StdError, ymax = Estimate + StdError), width = 0.1) +
  scale_color_manual(values = c("Non-significant" = "grey", "Significant" = "black")) +
  theme_classic() +
  theme(legend.position = 'none') +
  geom_hline(yintercept = 0, lty = 'dashed', col = 'grey') +
  labs(title = "Zero Hurdle Model Coefficient",
       y = "Estimate",
       x = "",
       color = "Significance")

# Display the plots


ggarrange(p_count, p_zero_hurdle, nrow = 2)


# 6. management effect: compare manag vs unmanaged stem densities, richness---------------------

create_plot <- function(df, x_var, y_var, y_label, colors) {
  df <- df %>% 
    dplyr::filter(cluster %in% cluster_id_keep) 
  
  ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) + 
    stat_summary(fun = "median", geom = "bar") +
    stat_summary(fun.data = function(z) {
      list(y = median(z), ymin = quantile(z, 0.25), ymax = quantile(z, 0.75))
    }, geom = "errorbar", width = .2) +
    ylab(y_label) +
    xlab('') +
    scale_fill_manual(values = colors) +
    theme_classic() + 
    theme(legend.position = 'none')
}

# Define colors for 'manag' categories
# Assuming 'manag' has two distinct categories
# Replace with actual category names and desired colors
colors <- c("Managed" = "darkolivegreen", "Unmanaged" = "darkolivegreen1")

# Using the function
p.richness <- create_plot(df_richness, "manag", "richness", "Richness", colors)
p.density <- create_plot(df_stems, "manag", "sum_stems", "Stems [ha]", colors)
p.vert <- create_plot(df_vert, "manag", "n_layers", "# Vertical Layers", colors)
p.IVI.max <- create_plot(df_IVI_max, "manag", "rIVI", "Dominance (rIVI [%])", colors)

# Arrange plots
ggarrange(p.density, p.vert, p.IVI.max, p.richness, ncol = 4)


# test significance - for 2 groups only
wilcox.test(sum_stems ~ manag, data = df_stems)
# test significance - for 2 groups only
wilcox.test(n_layers ~ manag, data = df_vert)

# test significance - for 2 groups only
wilcox.test(rIVI ~ manag, data = df_IVI_max)

# test significance - for 2 groups only
wilcox.test(richness ~ manag, data = df_richness)




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

          