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

library(cluster)   # cluster analysis






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
  filter(VegType == 'Saplings') %>%  # 20 cm - 2 m
  group_by(cluster) %>% 
  summarise(sum_stems_sapling = sum(total_stems))# %>% 

df_stems_mature <- stem_dens_ha_cluster_sum %>% 
  filter(VegType == 'Mature') %>% 
  group_by(cluster) %>% 
  summarise(sum_stems_mature = sum(total_stems))# %>% 

df_stems_juveniles <- stem_dens_ha_cluster_sum %>% 
  filter(VegType == 'Juveniles') %>%  # > 2m heigh
  group_by(cluster) %>% 
  summarise(sum_stems_juvenile = sum(total_stems))# %>% 



# get saplings and juveniles
df_stems_sapl_juveniles <- stem_dens_ha_cluster_sum %>% 
  filter(VegType != 'Mature') %>% 
  group_by(cluster, management_intensity) %>% 
  summarise(sum_stem_sapl_juven = sum(total_stems))# %>% 


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
    median_spei = median(spei, na.rm = TRUE),
    sd_spei = sd(spei, na.rm = TRUE)
  )

# <<Step 3: Plot - skip, long to run --------------------------
# fig_spei_all <- ggplot(spei_summary, aes(x = date, y = median_spei)) +
#   # Blue fill for positive values
#   geom_ribbon(data = subset(spei_summary, median_spei > 0),
#               aes(ymin = 0, ymax = median_spei),
#               fill = "blue", alpha = 1) +
#   # Red fill for negative values
#   geom_ribbon(data = subset(spei_summary, median_spei <= 0),
#               aes(ymin = median_spei, ymax = 0),
#               fill = "red", alpha = 1) +
#   theme_bw() +
#   labs(title = "SPEI median values by Scale",
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
  summarise(spei1 = median(spei1, na.rm = T),
            spei3 = median(spei3, na.rm = T),
            spei6 = median(spei6, na.rm = T),
            spei12 = median(spei12, na.rm = T),
            spei24 = median(spei24, na.rm = T)) #%>% 


spei_subplot_drought <- spei_wide %>% 
  ungroup(.) %>% 
  dplyr::filter(year %in% 2018 ) %>% # & month %in% 4:9# select just vegetation season
  dplyr::filter_all(all_vars(!is.infinite(.))) %>% # remove all infinite values
  #View()
  group_by(ID) %>%
  summarise(drought_spei1 = median(spei1, na.rm = T),
            drought_spei3 = median(spei3, na.rm = T),
            drought_spei6 = median(spei6, na.rm = T),
            drought_spei12 = median(spei12, na.rm = T),
            drought_spei24 = median(spei24, na.rm = T)) #%>% 

climate_subplot_drought <- climate %>% 
  ungroup(.) %>% 
  dplyr::filter(year %in% 2018 ) %>% # & month %in% 4:9# select just vegetation season
  dplyr::filter_all(all_vars(!is.infinite(.))) %>% # remove all infinite values
  #View()
  group_by(ID) %>%
  summarise(drought_tmp = median(tmp, na.rm = T),
            drought_prcp = median(prec, na.rm = T)) 


# merge all predictors together, ID (=subplot) level per year: having a median per years 2018-2023
df_predictors <- 
  climate %>% 
  left_join(climate_subplot_drought, by = join_by(ID)) %>% 
  left_join(dplyr::select(distance_to_edge, c(-country, -dist_ID ))) %>% 
  left_join(dplyr::select(disturbance_chars, c(-region, -country ))) %>% 
  left_join(soil) %>%
  left_join(spei_subplot) %>%          # clim values are medians from 2018-2023
  left_join(spei_subplot_drought) %>%  # clim values are medians from 2018-2020
  left_join(dplyr::select(terrain, c(-country, -region, -cluster.x, -cluster.y))) #%>% 
  #mutate(cluster = str_sub(ID, 4, -3)) 









length(unique(df_predictors$cluster)) # 957
anyNA(df_predictors)
# merge predcitors on cluster level: calculate medians
# keep only temp and prcp: 2021-2023 average
df_predictors_plot <- 
  df_predictors %>% 
  #dplyr::filter(year %in% 2018:2023) %>% 
  group_by(cluster) %>% 
  summarise(tmp          = median(tmp, na.rm = T),
            prec         = median(prec, na.rm = T),
            tmp_z        = median(tmp_z, na.rm = T),
            prcp_z       = median(prcp_z, na.rm = T),
            # Summary for drought years (2018-2020)
            drought_tmp  = median(drought_tmp, na.rm = TRUE),
            drought_prcp = median(drought_prcp , na.rm = TRUE),
            
            spei1        = median(spei1, na.rm = T),
            spei3        = median(spei3, na.rm = T),
            spei6        = median(spei6, na.rm = T),
            spei12       = median(spei12, na.rm = T),
            spei24       = median(spei24, na.rm = T),
            drought_spei1   = median(drought_spei1, na.rm = T),
            drought_spei3   = median(drought_spei3, na.rm = T),
            drought_spei6   = median(drought_spei6, na.rm = T),
            drought_spei12   = median(drought_spei12, na.rm = T),
            drought_spei24   = median(drought_spei24, na.rm = T),
            distance         = median(distance, na.rm = T),
            disturbance_year = as.integer(median(disturbance_year, na.rm = T)),
            disturbance_severity= median(disturbance_severity, na.rm = T),
            disturbance_agent= as.integer(median(disturbance_agent, na.rm = T)),
            elevation        = median(elevation, na.rm = T),
            sand_extract     = median(sand_extract, na.rm = T),
            silt_extract     = median(silt_extract, na.rm = T),
            clay_extract     = median(clay_extract, na.rm = T),
            depth_extract    = median(depth_extract, na.rm = T),
            av.nitro         = median(av.nitro, na.rm = T),
            slope            = median(slope, na.rm = T),
            aspect           = median(aspect, na.rm = T)
            )


 df_predictors_subplot <- 
  df_predictors %>% 
 # dplyr::filter(year %in% 2018:2023) %>% 
  group_by(ID) %>% 
  summarise(tmp = median(tmp, na.rm = T),
            prec = median(prec, na.rm = T),
            # Calculate drought temperature and precipitation only for 2018-2020
            drought_tmp = median(drought_tmp, na.rm = TRUE),
            drought_prcp = median(drought_prcp , na.rm = TRUE),
            tmp_z = median(tmp_z, na.rm = T),
            prcp_z = median(prcp_z, na.rm = T),
            spei1   = median(spei1, na.rm = T),
            spei3   = median(spei3, na.rm = T),
            spei6   = median(spei6, na.rm = T),
            spei12   = median(spei12, na.rm = T),
            spei24   = median(spei24, na.rm = T),
            distance = median(distance, na.rm = T),
            disturbance_year= as.integer(median(disturbance_year, na.rm = T)),
            disturbance_severity= median(disturbance_severity, na.rm = T),
            disturbance_agent= as.integer(median(disturbance_agent, na.rm = T)),
            elevation = median(elevation, na.rm = T),
            #region = median(tmp, na.rm = T),
            #country = median(tmp, na.rm = T),
            sand_extract= median(sand_extract, na.rm = T),
            silt_extract = median(silt_extract, na.rm = T),
            clay_extract = median(clay_extract, na.rm = T),
            depth_extract = median(depth_extract, na.rm = T),
            av.nitro    = median(av.nitro, na.rm = T),
            slope   = median(slope, na.rm = T),
            aspect= median(aspect, na.rm = T))


#fwrite(df_predictors_subplot, 'outData/all_predictors_subplot.csv')

fwrite(df_predictors_plot, 'outData/all_predictors_plot.csv')

# Get climate plots for map: TEMP, PREC, SPEI  -------------

# Custom function to return median and IQR for stat_summary
median_iqr <- function(y) {
  median_val <- median(y, na.rm = TRUE)
  iqr_val <- IQR(y, na.rm = TRUE)
  ymin <- median_val - (iqr_val / 2)
  ymax <- median_val + (iqr_val / 2)
  return(c(y = median_val, ymin = ymin, ymax = ymax))
}


# climate full
reference_period <- 1980:2010

# read temp and prec over years
climate_full           <- fread("outData/climate_1980_2023.csv")

spei_overview <- spei %>% 
  dplyr:: filter(!str_detect(ID, "^17")) %>% # remove Italy
  filter(!if_any(everything(), is.infinite)) %>% 
  group_by(cluster, year) %>% 
  dplyr::summarise(spei = median(spei, na.rm =T)) %>% 
  mutate(class = case_when(year %in% 2018:2020 ~ 'drought',
                           TRUE ~ 'ref')) #%>% 

# get simplified df
clim_overview <- climate_full %>% 
  dplyr:: filter(!str_detect(ID, "^17")) %>% # remove Italy
  mutate(cluster = str_sub(ID, 1, -3)) %>% 
  ungroup(.) %>% 
  group_by(cluster, year) %>% 
  dplyr::summarize(tmp = median(tmp),
                   prec = median(prec)) %>%
  mutate(class = case_when(year %in% 2018:2020 ~ 'drought',
                           TRUE ~ 'ref')) #%>% 

# get reference values
ref_spei <- median(spei_overview$spei[spei_overview$year %in% reference_period], na.rm =T)
ref_tmp  <- median(clim_overview$tmp[clim_overview$year %in% reference_period])
ref_prec <- median(clim_overview$prec[clim_overview$year %in% reference_period])


p.map.spei <- spei_overview %>% 
  ggplot(aes(x = year, y = spei, color = class)) +
  stat_summary(fun.data = median_iqr, geom = "pointrange", size = 0.1) +
  stat_summary(fun = median, geom = "point", size = 0.2) +
  geom_hline(yintercept = ref_spei, lty = 'dashed', col = 'grey70') +
  scale_color_manual(values = c('red', 'black')) +
  labs(x = "", y = expression(paste("SPEI [dim.]", sep=""))) +
  theme_classic() +
  theme(legend.position = "NULL")



# Plot for temperature
p.map.temp <- clim_overview %>% 
  ggplot(aes(x = year,
             y = tmp,
             color = class)) +
  stat_summary(fun.data = median_iqr, geom = "pointrange", size = 0.1) +  # Use median and IQR for pointrange
  stat_summary(fun = median, geom = "point", size = 0.2) +     # Add point for median
  geom_hline(yintercept = ref_tmp, lty = 'dashed', col = 'grey70') +
  scale_color_manual(values = c('red', 'black')) +
  labs(x = "", y =  expression(paste("Temperature [", degree, "C]", sep=""))) +
  theme_classic() +
  theme(legend.position = "NULL") 


# Plot for precipitation
p.map.prec <- clim_overview %>% 
  ggplot(aes(x = year,
             y = prec,
             color = class)) +
  stat_summary(fun.data = median_iqr, geom = "pointrange", size = 0.1) +  # Use median and IQR for pointrange
  stat_summary(fun = median, geom = "point", size = 0.2) +     # Add point for median
  geom_hline(yintercept = ref_prec, lty = 'dashed', col = 'grey70') +
  scale_color_manual(values = c('red', 'black')) +
  labs(x = "", y =  expression(paste("Precipitation [mm]", sep=""))) +
  theme_classic() +
  theme(legend.position = "NULL")

windows(7,2)
p.clim.map <- ggarrange(p.map.temp, p.map.prec,p.map.spei, ncol = 3)
p.clim.map
ggsave(filename = 'outFigs/Fig1.png', plot = p.clim.map, width = 7, 
       height = 2, dpi = 300, bg = 'white')




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

#length(unique(df_plot$cluster))
length(unique(df_plot_full$cluster)) # 849!   - final clusters, 4-5 plots
length(unique(df_predictors_plot$cluster))  # 957 - all clusters, from even with less plots




# Environmental cluster analysis --------------------------------------------------
# inspect if drought in 2018 would not be a better indicator to get clim clusters

df_fin <- df_fin %>% 
  rename(prcp = prec) %>% 
  rename(site = cluster)

## Cluster: Climate-environment: SPEI 3 -----------------------------------------------
# for spei3 - mean per 2018-2023 

# Subset the relevant columns
data_subset_clim <- df_fin[, c("tmp", "prcp", "tmp_z", "prcp_z", "spei3", "sand_extract", "clay_extract", "depth_extract", "av.nitro")]
#
# Standardize the data
data_scaled_clim <- scale(data_subset_clim)

# Find which number of clusters is teh best
# Perform K-means clustering for different values of k
set.seed(3)
max_clusters <- 10
sil_width <- numeric(max_clusters)
for (k in 2:max_clusters) {
  kmeans_result_clim <- kmeans(data_scaled_clim, centers = k, nstart = 25)
  sil <- silhouette(kmeans_result_clim$cluster, dist(data_scaled_clim))
  sil_width[k] <- mean(sil[, 3])
}

# Plot Silhouette width for different values of k
plot(1:max_clusters, sil_width, type = "b", xlab = "Number of clusters", ylab = "Average Silhouette width", main = "Silhouette Analysis for K-means Clustering")

# Determine the optimal number of clusters
optimal_k_clim <- 3  # from Kilian's study

# Perform K-means clustering with the optimal number of clusters
set.seed(3)
kmeans_result <- kmeans(data_scaled_clim, centers = optimal_k_clim, nstart = 25)

# Add cluster assignments to the original data
df_fin$clim_cluster_spei3 <- kmeans_result$cluster

# Perform PCA for visualization - use Pc1 and Pc1
pca_result_clim <- prcomp(data_scaled_clim)

# plot PCA results with the most important variables: 
biplot(pca_result_clim, main = "PCA Biplot")

#pairs(data_scaled_clim)

# Convert to data frame if it's not already
data_scaled_clim_df <- as.data.frame(data_scaled_clim)

# Create a ggpairs plot
#ggpairs(data_scaled_clim_df)

# drought 


# classify based on scatter plots and groups: 

# 1 - cold, wet - similar soil coonditions to cluster 3
# 2 - hot, dry, sandy - more sand, less clay then 3, more av.nitro
# 3 - drier, mild, clay - more clay tehn n3

df_fin <- df_fin %>% 
  mutate(clim_class = case_when(
    clim_cluster_spei3 == 1 ~ "wet-warm-clay",  # Cluster 1: wet, cold, clay
    clim_cluster_spei3 == 2 ~ "hot-dry-sand",  # Cluster 2: hot, dry, sandy (more sand, less clay, more av.nitro than cluster 3)
    clim_cluster_spei3 == 3 ~ "hot-dry-clay"    # Cluster 3: hot, dry, more clay
  )) %>% 
  mutate( clim_cluster_spei3 = as.factor(clim_cluster_spei3),
          clim_class = as.factor(clim_class))




# Make a scatetr plot of groups: ------------------------
# scatter plot: tmp vs spei

# Calculate the mean and standard deviation for each cluster
median_iqr_df <- df_fin %>%
  group_by(clim_class) %>%
  summarise(
    median_tmp = median(tmp, na.rm = TRUE),
    iqr_tmp = IQR(tmp, na.rm = TRUE),
    median_spei = median(spei3, na.rm = TRUE),
    iqr_spei = IQR(spei3, na.rm = TRUE)
  )



# Create the scatter plot with ellipses, mean points, and error bars
fig_spei_tmp_clusters <- ggplot(data = df_fin, aes(x = tmp, 
                          y = spei3,
                          group = clim_class,
                          color = clim_class, 
                          fill = clim_class#,
                          # size = prcp
)) +  # Map point size to prcp
  geom_point(aes(size = prcp), shape = 16, alpha = 0.5) +  # No fixed size here
  scale_size_continuous(
    range = c(0.01, 1.7),  # Adjust the size range of points
    name = "Precipitation [mm]",
    breaks =  c(550, 1700) #,
      #c(round(min(df_fin$prcp, na.rm = TRUE),0), 
      #         round(max(df_fin$prcp, na.rm = TRUE),0))
    ) + # Add a size legend for precipitation
  stat_ellipse(aes(group = clim_class), type = "norm", alpha = 0.3, geom = "polygon") +
  geom_point(data = median_iqr_df, aes(x = median_tmp, y = median_spei), size = 2.8, shape = 16, color = 'white') +
  geom_point(data = median_iqr_df, aes(x = median_tmp, y = median_spei), size = 1.5, shape = 16) +
  geom_errorbar(data = median_iqr_df, aes(x = median_tmp, ymin = median_spei - iqr_spei/2, 
                                          ymax = median_spei + iqr_spei/2, y = median_spei), 
                width = 0.2) +
  geom_errorbarh(data = median_iqr_df, aes(y = median_spei, xmin = median_tmp - iqr_tmp/2, xmax = median_tmp + iqr_tmp/2, x = median_tmp), 
                 height = 0.01) +
  theme_classic2() +
  labs(title = "",
       x = expression(Temperature ~ "[" * degree * C * "]"),
       y = "SPEI3 [dim.]",
       fill = 'Clim-ENV cluster',
       color = 'Clim-ENV cluster') +
  scale_fill_manual(values = c("orange", "red", "blue")) +  # Customize fill colors
  scale_color_manual(values = c("orange", "red", "blue")) +  # Customize point colors
  # Set all text sizes to 8 using theme()
  theme(
    text = element_text(size = 8),             # Set the base text size to 8
    axis.text = element_text(size = 8),        # Set axis text size to 8
    axis.title = element_text(size = 8),       # Set axis title text size to 8
    legend.text = element_text(size = 8),      # Set legend text size to 8
    legend.title = element_text(size = 8),     # Set legend title text size to 8
    plot.title = element_text(size = 8)        # Set plot title text size to 8
  )


(fig_spei_tmp_clusters )
ggsave(filename = 'outFigs/fig_clim_clusters.png', plot = fig_spei_tmp_clusters, 
       width = 4.5, height = 3, dpi = 300, 
       bg = 'white')


# GEt summary table of another variables in 3 env clustersL

# Summarize the data based on 'clim_cluster_spei3' for the specified variables
# Summarize the data based on 'clim_cluster_spei3' for the specified variables
summary_table <- df_fin %>%
  group_by(clim_class) %>%
  summarise(
    temperature = paste0(round(median(tmp, na.rm = TRUE), 2), " (IQR: ", round(IQR(tmp, na.rm = TRUE), 2), ")"),
    precipitation = paste0(round(median(prcp, na.rm = TRUE), 2), " (IQR: ", round(IQR(prcp, na.rm = TRUE), 2), ")"),
    spei = paste0(round(median(spei3, na.rm = TRUE), 2), " (IQR: ", round(IQR(spei3, na.rm = TRUE), 2), ")"),
    sand = paste0(round(median(sand_extract, na.rm = TRUE), 2), " (IQR: ", round(IQR(sand_extract, na.rm = TRUE), 2), ")"),
    clay = paste0(round(median(clay_extract, na.rm = TRUE), 2), " (IQR: ", round(IQR(clay_extract, na.rm = TRUE), 2), ")"),
    depth = paste0(round(median(depth_extract, na.rm = TRUE), 2), " (IQR: ", round(IQR(depth_extract, na.rm = TRUE), 2), ")"),
    av_nitrogen = paste0(round(median(av.nitro, na.rm = TRUE), 2), " (IQR: ", round(IQR(av.nitro, na.rm = TRUE), 2), ")")
  )

# Print the summary table
print(summary_table)


sjPlot::tab_df(summary_table,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = F,
               file="outTable/clim_cluster_summary_env_conditions.doc",
               digits = 1) 
