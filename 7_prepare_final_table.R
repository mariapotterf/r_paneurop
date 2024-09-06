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
library(tidyr)
library(terra)
library(ggplot2)
library(lubridate)
library(stringr)

library(ggpubr)

#library(dunn.test) # Duncan test, post hoc Kruskal-Wallis






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

# create idividual table for each of the vertical regeneration class 
df_stems_saplings <- stem_dens_ha_cluster_sum %>% 
  filter(VegType == 'Saplings') %>% 
  group_by(cluster) %>% 
  summarise(sum_stems_sapling = sum(total_stems))# %>% 

df_stems_mature <- stem_dens_ha_cluster_sum %>% 
  filter(VegType == 'Mature') %>% 
  group_by(cluster) %>% 
  summarise(sum_stems_mature = sum(total_stems))# %>% 

df_stems_juveniles <- stem_dens_ha_cluster_sum %>% 
  filter(VegType == 'Juveniles') %>% 
  group_by(cluster) %>% 
  summarise(sum_stems_juvenile = sum(total_stems))# %>% 



# get saplings and juveniles
df_stems_sapl_juveniles <- stem_dens_ha_cluster_sum %>% 
  filter(VegType != 'Mature') %>% 
  group_by(cluster, management_intensity) %>% 
  summarise(sum_stem_sapl_juven = sum(total_stems))# %>% 


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


# SPEI - development over years by scale (number of months)
# Step 1: Create a Date Column
spei_date <- spei %>%
  mutate(date = ymd(paste(year, month, "01", sep = "-")))  # Create a date column

# Step 2: Calculate Summary Statistics
spei_summary <- spei_date %>%
  group_by(scale, date) %>%
  summarise(
    mean_spei = mean(spei, na.rm = TRUE),
    sd_spei = sd(spei, na.rm = TRUE)
  )

# <<Step 3: Plot - skip, long to run --------------------------
# fig_spei_all <- ggplot(spei_summary, aes(x = date, y = mean_spei)) +
#   # Blue fill for positive values
#   geom_ribbon(data = subset(spei_summary, mean_spei > 0),
#               aes(ymin = 0, ymax = mean_spei),
#               fill = "blue", alpha = 1) +
#   # Red fill for negative values
#   geom_ribbon(data = subset(spei_summary, mean_spei <= 0),
#               aes(ymin = mean_spei, ymax = 0),
#               fill = "red", alpha = 1) +
#   theme_bw() +
#   labs(title = "SPEI Mean values by Scale",
#        x = "Year",
#        y = "SPEI Value",
#        color = "Scale",
#        fill = "Scale") +
#   theme(legend.position = "bottom") +
#   facet_grid(scale~.)
# 


# ggsave(filename = 'outFigs/fig_spei_all.png', 
#        plot = fig_spei_all, 
#        width = 7, height = 7, dpi = 300, bg = 'white')
# >>>>



# Convert to wide format using pivot_wider()
spei_wide <- spei %>%
  mutate(scale = paste0("spei", scale)) %>%  # Create column names like spei1, spei3, etc.
  pivot_wider(names_from = scale, values_from = spei)




spei_subplot <- spei_wide %>% 
  ungroup(.) %>% 
  dplyr::filter(year %in% 2018:2023 ) %>% # & month %in% 4:9# select just vegetation season
  dplyr::filter_all(all_vars(!is.infinite(.))) %>% # remove all infinite values
  #View()
  group_by(ID, year) %>%
  summarise(spei1 = mean(spei1, na.rm = T),
            spei3 = mean(spei3, na.rm = T),
            spei6 = mean(spei6, na.rm = T),
            spei12 = mean(spei12, na.rm = T),
            spei24 = mean(spei24, na.rm = T)) #%>% 


spei_subplot_drought <- spei_wide %>% 
  ungroup(.) %>% 
  dplyr::filter(year %in% 2018:2020 ) %>% # & month %in% 4:9# select just vegetation season
  dplyr::filter_all(all_vars(!is.infinite(.))) %>% # remove all infinite values
  #View()
  group_by(ID, year) %>%
  summarise(drought_spei1 = mean(spei1, na.rm = T),
            drought_spei3 = mean(spei3, na.rm = T),
            drought_spei6 = mean(spei6, na.rm = T),
            drought_spei12 = mean(spei12, na.rm = T),
            drought_spei24 = mean(spei24, na.rm = T)) #%>% 




spei_plot <- spei %>% 
  ungroup(.) %>% 
  dplyr::filter(year %in% 2018:2023 ) %>% # & month %in% 4:9# select just vegetation season
  dplyr::filter_all(all_vars(!is.infinite(.))) %>% # remove all infinite values
  #View()
  group_by(cluster, year) %>%
  summarise(spei = mean(spei, na.rm = T)) #%>% 



# merge all predictors together, ID (=subplot) level
df_predictors <- 
  climate %>% 
  left_join(dplyr::select(distance_to_edge, c(-country, -dist_ID ))) %>% 
  left_join(dplyr::select(disturbance_chars, c(-region, -country ))) %>% 
  left_join(soil) %>%
  left_join(spei_subplot) %>%          # clim values are means from 2018-2023
  left_join(spei_subplot_drought) %>%  # clim values are means from 2018-2020
  
  left_join(dplyr::select(terrain, c(-country, -region, -cluster.x, -cluster.y))) #%>% 
  #mutate(cluster = str_sub(ID, 4, -3)) 

length(unique(df_predictors$cluster)) # 957
anyNA(df_predictors)
# merge predcitors on cluster level: calculate means
# keep only temp and prcp: 2021-2023 average
df_predictors_plot <- 
  df_predictors %>% 
  #dplyr::filter(year %in% 2018:2023) %>% 
  group_by(cluster) %>% 
  summarise(tmp = mean(tmp, na.rm = T),
            prec = mean(prec, na.rm = T),
            tmp_z = mean(tmp_z, na.rm = T),
            prcp_z = mean(prcp_z, na.rm = T),
            spei1   = mean(spei1, na.rm = T),
            spei3   = mean(spei3, na.rm = T),
            spei6   = mean(spei6, na.rm = T),
            spei12   = mean(spei12, na.rm = T),
            spei24   = mean(spei24, na.rm = T),
            drought_spei1   = mean(drought_spei1, na.rm = T),
            drought_spei3   = mean(drought_spei3, na.rm = T),
            drought_spei6   = mean(drought_spei6, na.rm = T),
            drought_spei12   = mean(drought_spei12, na.rm = T),
            drought_spei24   = mean(drought_spei24, na.rm = T),
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


 df_predictors_subplot <- 
  df_predictors %>% 
 # dplyr::filter(year %in% 2018:2023) %>% 
  group_by(ID) %>% 
  summarise(tmp = mean(tmp, na.rm = T),
            prec = mean(prec, na.rm = T),
            tmp_z = mean(tmp_z, na.rm = T),
            prcp_z = mean(prcp_z, na.rm = T),
            spei1   = mean(spei1, na.rm = T),
            spei3   = mean(spei3, na.rm = T),
            spei6   = mean(spei6, na.rm = T),
            spei12   = mean(spei12, na.rm = T),
            spei24   = mean(spei24, na.rm = T),
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


#fwrite(df_predictors_subplot, 'outData/all_predictors_subplot.csv')

fwrite(df_predictors_plot, 'outData/all_predictors_plot.csv')

# Get climate plots for map: TEMP, PREC, SPEI  -------------

# climate full
reference_period <- 1980:2010

# read temp and prec over years
climate_full           <- fread("outData/climate_1980_2023.csv")

spei_overview <- spei %>% 
  dplyr:: filter(!str_detect(ID, "^17")) %>% # remove Italy
  filter(!if_any(everything(), is.infinite)) %>% 
  group_by(cluster, year) %>% 
  dplyr::summarise(spei = mean(spei, na.rm =T)) %>% 
  mutate(class = case_when(year %in% 2018:2020 ~ 'drought',
                           TRUE ~ 'ref')) #%>% 

ref_spei <- mean(spei_overview$spei[spei_overview$year %in% reference_period], na.rm =T)

p.map.spei <- spei_overview %>% 
  ggplot(aes(x = year,
             y = spei,
             color = class)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange") +
  # This will add a point for the mean of Local Variance
  stat_summary(fun = mean, geom = "point", size = 0.7) +
  geom_hline(yintercept = ref_spei, lty = 'dashed', col = 'grey70') +
  scale_color_manual(values = c('red', 'black')) +
  labs(x = "", y =  expression(paste("SPEI [dim.]", sep=""))) +
  theme_classic() +
  theme(legend.position = "NULL") 


# get simplified df
clim_overview <- climate_full %>% 
  dplyr:: filter(!str_detect(ID, "^17")) %>% # remove Italy
  mutate(cluster = str_sub(ID, 1, -3)) %>% 
  ungroup(.) %>% 
  group_by(cluster, year) %>% 
  dplyr::summarize(tmp = mean(tmp),
                   prec = mean(prec)) %>%
  mutate(class = case_when(year %in% 2018:2020 ~ 'drought',
                           TRUE ~ 'ref')) #%>% 

# make plots for map:
ref_tmp <- mean(clim_overview$tmp[clim_overview$year %in% reference_period])
ref_prec <- mean(clim_overview$prec[clim_overview$year %in% reference_period])

p.map.temp <- clim_overview %>% 
  ggplot(aes(x = year,
             y = tmp,
             color = class)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange") +
  # This will add a point for the mean of Local Variance
  stat_summary(fun = mean, geom = "point", size = 0.7) +
  geom_hline(yintercept = ref_tmp, lty = 'dashed', col = 'grey70') +
  scale_color_manual(values = c('red', 'black')) +
  labs(x = "", y =  expression(paste("Temperature [", degree, "C]", sep=""))) +
  theme_classic() +
  theme(legend.position = "NULL") 


p.map.prec <- clim_overview %>% 
  ggplot(aes(x = year,
             y = prec,
             color = class)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange") +
  # This will add a point for the mean of Local Variance
  stat_summary(fun = mean, geom = "point", size = 0.7) +
  geom_hline(yintercept = ref_prec , lty = 'dashed', col = 'grey70') +
  scale_color_manual(values = c('red', 'black')) +
  labs(x = "", y =  expression(paste("Precipitation [mm]", sep=""))) +
  theme_classic() +
  theme(legend.position = "NULL") 

windows(7,2)
p.clim.map <- ggarrange(p.map.temp, p.map.prec,p.map.spei, ncol = 3)
p.clim.map
ggsave(filename = 'outFigs/Fig1.png', plot = p.clim.map, width = 7, height = 2, dpi = 300, bg = 'white')




# Merge predictors with vegetation data ---------------------------------------------------
# select the dominant species per cluster
df_IVI_max <- df_IVI %>% 
  dplyr::select(cluster, Species, country,rIVI) %>% 
  ungroup(.) %>% 
  group_by(cluster, country) %>% 
  filter(rIVI == max(rIVI)) %>% 
  slice(1)

# plot level
df_plot_veg <-
  df_IVI_max %>% 
  left_join(df_richness) %>%
  left_join(df_vert) %>%
  left_join(df_stems) %>%
  left_join(df_stems_juveniles) %>%
  left_join(df_stems_saplings ) %>% 
  left_join(df_stems_mature)

anyNA(df_plot_veg)


# plot data share with predictors:
# rename not necessary predictors for climate clustering, rename veg variables
df_plot_full <- df_plot_veg %>% 
  left_join(df_predictors_plot, by = join_by(cluster)) %>% 
  ungroup() %>% 
  dplyr::select(-c(disturbance_year, 
                   disturbance_agent,
                   silt_extract, #country,
                   elevation,
                   slope, aspect)) %>% 
  dplyr::rename(dominant_species = Species, 
                stem_density     = sum_stems,
                distance_edge    = distance,
                n_vertical       = n_layers) %>% 
  left_join(dat_manag_intensity_cl, by = join_by(cluster, management_intensity))# %>% 
  
# set up proper structure (df) fr analysis  and claim characters as factors
df_fin <- as.data.frame(df_plot_full) %>% 
  mutate(cluster = factor(cluster),
         dominant_species = factor(dominant_species),
         country = factor(country))


plot_n <- length(unique(df_plot_full$cluster))


fwrite(df_fin, 'outData/indicators_for_cluster_analysis.csv')

#length(unique(df_plot$cluster))
length(unique(df_plot_full$cluster)) # 849!   - final clusters, 4-5 plots
length(unique(df_predictors_plot$cluster))  # 957 - all clusters, from even with less plots

