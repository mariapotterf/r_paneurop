# create final database
# merge veg data
# predictors - climate, disturbance chars, soils
# dependent - regeneretion, density, richness
# output file
# 

gc()

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


### read data --------------------------------------------------------------------

climate           <- fread("outData/climate_18_23.csv")
spei              <- fread("outData/SPEI.csv")
distance_to_edge  <- fread("outData/distance_to_edge.csv")
disturbance_chars <- fread("outData/disturbance_chars.csv")
soil              <- as.data.frame(terra::vect("rawData/extracted_soil_data/extracted_soil_data_completed.gpkg"))
terrain           <- fread("outData/terrain.csv")

# get vegetation data
load("outData/veg.Rdata")


### clean up data -----------------------------------------------------------------
# sum stems first up across vertical groups!!

# get only samplings
df_stems_sapl <- stem_dens_ha_cluster_sum %>% 
  filter(VegType == 'Saplings') %>% 
  group_by(cluster, management_intensity) %>% 
  summarise(sum_stems = sum(total_stems))# %>% 

# get saplings and juveniles
df_stems_both <- stem_dens_ha_cluster_sum %>% 
  filter(VegType != 'Mature') %>% 
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
                n_vertical = n_layers) %>% 
  left_join(dat_manag_intensity_cl, by = join_by(cluster, management_intensity))# %>% 
  
  

fwrite(df_cluster_full, 'outData/indicators_for_cluster_analysis.csv')

length(unique(df_cluster$cluster))
length(unique(df_cluster_full$cluster)) # 849!   - final clusters, 4-5 plots
length(unique(df_predictors_cluster$cluster))  # 957 - all clusters, from even with less plots

cluster_n <- length(unique(df_cluster_full$cluster))


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


#View(df_cluster_full)
# get summary for the scatter plot, one point is one country & management
df_summary <- df_cluster_full %>%
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
p1 <- df_cluster_full %>% 
  ggplot(aes(x = management_intensity ,
             y = rIVI))+ 
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = 'loess') +
  theme_bw() +
  theme(aspect.ratio = 1) 
  
  
p2 <- df_cluster_full %>% 
    ggplot(aes(x = management_intensity,
               y = richness))+ 
  geom_jitter(alpha = 0.5) +
    geom_smooth(method = 'loess')+
  theme_bw() +
  theme(aspect.ratio = 1) 

p3 <- df_cluster_full %>% 
  ggplot(aes(x = management_intensity,
             y = stem_density))+ 
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = 'loess')+
  theme_bw() +
  theme(aspect.ratio = 1) 

p4 <- df_cluster_full %>% 
  ggplot(aes(x = management_intensity,
             y = n_vertical))+
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = 'loess')+
 
  theme_bw() +
  theme(aspect.ratio = 1) 
  


ggarrange(p1, p2, p3,p4, nrow = 2, ncol = 2, align = 'hv')


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




###### prepare plot with raster -------------------------------

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

##### Density plots: each variable against another one --------------------
# 1. rIVI vs stem_density
# 2. rIVI vs richness
# 3. rIVI vs n_vertical
# 4. stem_density vs richness
# 5. stem_density vs n_vertical
# 6. richness vs n_vertical 


######### rIVI vs sum_stems 
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     x = 'stem_density',
                                     y = 'rIVI')


p_IVI_stems <- make_2D_plot(df_plot_data, 
             x_lab = "Stem density [n/ha]",
             y_lab = "Importance value [%]")

(p_IVI_stems)
####### rIVI vs richnes
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     x = 'rIVI', 
                                     y = 'richness')


p_IVI_richness <- make_2D_plot(df_plot_data, 
                            x_lab = "Importance value [%]", # Default axis labels, can be customized
                            y_lab = "Richness [counts]"  )


####### rIVI vs vertical 
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     x = 'rIVI', 
                                     y = 'n_vertical')


p_IVI_vert <- make_2D_plot(df_plot_data, 
                               x_lab = "Importance value [%]", # Default axis labels, can be customized
                               y_lab = "# vertical layers [counts]"  )



####### stem_density vs richness 
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     x = 'stem_density', 
                                     y = 'richness')


p_dens_richness <- make_2D_plot(df_plot_data, 
                           x_lab = "Stem density [n/ha]", # Default axis labels, can be customized
                           y_lab = "Richness [counts]"  )



####### stem_density vs n_vertcal 
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     x = 'stem_density', 
                                     y = 'n_vertical')


p_dens_vertical <- make_2D_plot(df_plot_data, 
                                x_lab = "Stem density [n/ha]", # Default axis labels, can be customized
                                y_lab = "# vertical layers [counts]"  )


####### n_vertical vs richness 
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
          p_vertical_richness, nrow = 2, ncol = 3, common.legend = T, legend = 'bottom', align = 'hv',
          labels = c("[a]","[b]","[c]","[d]","[e]","[f]"))
#### DEnsity plot wth management 
###### rIVI vs management 
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     y = 'rIVI', 
                                     x = 'management_intensity')


p_IVI_manag <- make_2D_plot(df_plot_data, 
                            y_lab = "Importance value [%]", # Default axis labels, can be customized
                            x_lab = "Management intensity [%]"  )

(p_IVI_manag)

###### richness vs management 
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     y = 'richness', 
                                     x = 'management_intensity')


p_richness_manag <- make_2D_plot(df_plot_data, 
                            y_lab = "Richness [#]", # Default axis labels, can be customized
                            x_lab = "Management intensity [%]"  )

(p_richness_manag)



###### stem_density vs management
df_plot_data <- prepare_density_data(df_cluster_full, 
                                     y = 'stem_density', 
                                     x = 'management_intensity')


p_stem_dens_manag <- make_2D_plot(df_plot_data, 
                                 y_lab = "Stem density [n/ha]", # Default axis labels, can be customized
                                 x_lab = "Management intensity [%]"  )

(p_stem_dens_manag)




###### n_vertical vs management 
df_plot_data <- prepare_density_data(df_cluster_full, 
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




##### get summary table per quantiles -------------------------------------------------

# Represent results using quantiles, as they are skewed?
qntils = c(0, 0.01, 0.25, 0.5, 0.75, 0.90, 1)

qs_dat_fin <- 
  df_cluster_full %>% 
  ungroup(.) %>% 
  #filter(year %in% 2015:2021) %>% 
  dplyr::reframe(stem_density    = quantile(stem_density , qntils, na.rm = T ),
                 n_vertical      = quantile(n_vertical, qntils, na.rm = T ),
                 rIVI            = quantile(rIVI, qntils, na.rm = T ),
                 richness        = quantile(richness, qntils, na.rm = T )#,
                 ) %>% 
  t() %>%
  round(1) %>% 
  as.data.frame()

(qs_dat_fin)  


means_dat_fin <- 
  df_cluster_full %>% 
  ungroup(.) %>% 
  dplyr::reframe(stem_density    = mean(stem_density , na.rm = T ),
                 n_vertical      = mean(n_vertical, na.rm = T ),
                 rIVI            = mean(rIVI,na.rm = T ),
                 richness        = mean(richness, na.rm = T )
                 ) %>% 
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










# Descriptive tables  -----------------------------------------------------------

# Richness desc 

# decrsibe richness: how many clusters have how many species?
df_cluster_full %>% 
  group_by(richness) %>% 
  dplyr::summarize(count = n(),
                   prop = count/cluster_n*100)


# what are dominant species??
df_cluster_full %>% 
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
View(df_stock_density)

# make a table for potential regeneration delay - link with climate data
# make a barplots
df_regen_delay <- df_cluster_full %>% 
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
                        #x_var, 
                        y_var, 
                        #x_label = "", 
                        y_label = y_var,
                       x_annotate = 0, lab_annotate = "lab ann") {
  
  p <- ggplot(data, aes_string(x = "reg_delay_class", 
                          y = y_var)) +
    stat_summary(fun.data = median_iqr, geom = "pointrange") +
    stat_summary(fun = median, geom = "point") +
   # stat_summary(fun.data = "mean_se", geom = "pointrange") + # Specify the summary function explicitly
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






# Using the function to create plots
p1 <- create_plot(df_regen_delay, #x_var = "reg_delay_class", 
                  y_var = "tmp",  x_annotate = 1.5, lab_annotate = "0.06")
p2 <- create_plot(df_regen_delay, 
                  #x_var = "reg_delay_class", 
                  y_var = "prec",  x_annotate = 1.5, lab_annotate = "***")
p3 <- create_plot(df_regen_delay, 
                  #x_var = "reg_delay_class", 
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

df_cluster_full
windows()
pairs(stem_density    ~ tmp + tmp_z + prec + prcp_z + management_intensity+ salvage_intensity + protection_intensity + distance_edge, df_cluster_full)

pairs(stem_density    ~  management_intensity+ salvage_intensity + protection_intensity, df_cluster_full)


# test corelations 

plot(df_cluster_full$tmp, df_cluster_full$spei)
plot(df_cluster_full$prec, df_cluster_full$spei)
plot(df_cluster_full$tmp, df_cluster_full$prec)

plot( df_cluster_full$management_intensity, df_cluster_full$stem_density)
plot( df_cluster_full$salvage_intensity, df_cluster_full$stem_density)
plot( df_cluster_full$protection_intensity, df_cluster_full$stem_density)


# keep only spei instead of the tmp and prec

# keep also annomalies?
plot(df_cluster_full$tmp_z, df_cluster_full$spei)

# test, which one of teh variables exaplin teh stem density better?




cor(df_cluster_full$tmp_z, df_cluster_full$spei, method = 'spearman')
cor(df_cluster_full$tmp, df_cluster_full$spei, method = 'spearman')
#[1] -0.7187752
cor(df_cluster_full$prec, df_cluster_full$spei, method = 'spearman')
# [1] 0.2349638

# check how management is correlated?
cor(df_cluster_full$management_intensity, df_cluster_full$salvage_intensity, method = 'spearman')
cor(df_cluster_full$protection_intensity, df_cluster_full$salvage_intensity, method = 'spearman')



# decide: which parameters are better: temp_z, temp, prcp, prcp_z or spei?? or their interation?

# Calculate correlation matrix for predictors
correlation_matrix <- cor(df_cluster_full[, c('tmp_z', 'tmp', 'prec', 'prcp_z', 'spei')])


# remove spei
correlation_matrix <- cor(df_cluster_full[, c('tmp_z', 'tmp', 'prec', 'prcp_z')], method = 'spearman')

correlation_matrix




#5. test model ---------------------------------

df_model <- df_cluster_full %>% 
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
ggplot(df_model, aes(x = spei, y = stem_density)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
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

glmm1 <- glmmTMB(stem_density ~ tmp  + tmp_z + prcp_z + prec +
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


# use zero inflated model:
m.zi1 <- glmmTMB(stem_density ~ tmp  + tmp_z + prcp_z + prec + spei + tmp:spei +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   #country+
                   sand_extract + depth_extract + av.nitro +
                   (1 | country), 
                 ziformula = ~1,
                 data = df_model, 
                 family = nbinom2, na.action = "na.fail")


m.zi2 <- glmmTMB(stem_density ~ tmp  + tmp_z + prcp_z + 
                   prec + #remove, as it has lower effect then prc_z 
                   spei + tmp:spei +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   #country+
                   #sand_extract + depth_extract + av.nitro +
                   (1 | country), 
                 ziformula = ~1,
                 data = df_model, 
                 family = nbinom2, na.action = "na.fail")

# remove spei, get interaction between tmp_z and prcp_z
m.zi3 <- glmmTMB(stem_density ~ tmp  +# tmp_z + prcp_z + 
                   prec + #remove, as it has lower effect then prc_z 
                   #spei + 
                   tmp_z*prcp_z +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   #country+
                   #sand_extract + depth_extract + av.nitro +
                   (1 | country), 
                 ziformula = ~1,
                 data = df_model, 
                 family = nbinom2, na.action = "na.fail")

# Chek for multicollinearity
library(performance)

# Assuming your model is an 'lme4' or 'glmmTMB' object
(vif_values <- performance::check_collinearity(m.zi7))

# remove highly correlated values: tmp_z*prcp_z
m.zi4 <- glmmTMB(stem_density ~ #tmp*prec  +# tmp_z + prcp_z + 
                   #prec + #remove, as it has lower effect then prc_z 
                   #spei + 
                   tmp_z*prcp_z +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   #country+
                   #sand_extract + depth_extract + av.nitro +
                   (1 | country), 
                 ziformula = ~1,
                 data = df_model, 
                 family = nbinom2, na.action = "na.fail")


m.zi5 <- glmmTMB(stem_density ~ tmp  + tmp_z + prcp_z + 
                   prec + #remove, as it has lower effect then prc_z 
                   #spei + 
                   tmp:prec +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   #country+
                   #sand_extract + depth_extract + av.nitro +
                   (1 | country), 
                 ziformula = ~1,
                 data = df_model, 
                 family = nbinom2, na.action = "na.fail")


m.zi6 <- glmmTMB(stem_density ~ #tmp  + 
                   tmp_z*prcp_z + 
                   #prec + #remove, as it has lower effect then prc_z 
                   #spei + 
                   #tmp:prec +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   #country+
                   #sand_extract + depth_extract + av.nitro +
                   (1 | country), 
                 ziformula = ~1,
                 data = df_model, 
                 family = nbinom2, na.action = "na.fail")

# remove high VIF values
m.zi7 <- glmmTMB(stem_density ~ #tmp  + 
                   #tmp_z + 
                   #prec + #remove, as it has lower effect then prc_z 
                   spei + 
                   #tmp:prec +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   #country+
                   #sand_extract + depth_extract + av.nitro +
                   (1 | country), 
                 ziformula = ~1,
                 data = df_model, 
                 family = nbinom2, na.action = "na.fail")

# keep only one of climatic predictors
m.zi8 <- glmmTMB(stem_density ~ #tmp  + 
                   tmp_z + 
                   #prec + #remove, as it has lower effect then prc_z 
                   #spei + 
                   #tmp:prec +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   #country+
                   #sand_extract + depth_extract + av.nitro +
                   (1 | country), 
                 ziformula = ~1,
                 data = df_model, 
                 family = nbinom2, na.action = "na.fail")


m.zi9 <- glmmTMB(stem_density ~ prcp_z + #tmp  + 
                   #tmp_z + 
                   #prec + #remove, as it has lower effect then prc_z 
                   #spei + 
                   #tmp:prec +
                   disturbance_severity + distance_edge  + 
                   management_intensity + 
                   clay_extract + 
                   #country+
                   #sand_extract + depth_extract + av.nitro +
                   (1 | country), 
                 ziformula = ~1,
                 data = df_model, 
                 family = nbinom2, na.action = "na.fail")



AIC(m.zi1, m.zi2,m.zi3,m.zi4,m.zi5,m.zi6,m.zi7,m.zi8,m.zi9)
summary(m.zi9)
plot(allEffects(m.zi9))

simulateResiduals(m.zi9, plot = T)

AIC(global.negbin4, glmm1,glmm2,glmm3,glmm4)
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

          