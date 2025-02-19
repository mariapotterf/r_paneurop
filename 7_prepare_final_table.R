# create final database
# merge veg data
# predictors - climate, disturbance chars, soils - update climate predictors for Germany from DWD data
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
library(cluster)   # cluster analysis






# Read data --------------------------------------------------------------------

# all states - from ERA, 10 km res, level of subplots
climate           <- fread("outData/climate_18_23.csv")
spei              <- fread("outData/SPEI.csv")
distance_to_edge  <- fread("outData/distance_to_edge.csv")
disturbance_chars <- fread("outData/disturbance_chars.csv")
soil              <- as.data.frame(terra::vect("rawData/extracted_soil_data/extracted_soil_data_completed.gpkg"))
terrain           <- fread("outData/terrain.csv")

# higher resolution climate data: germany (1 km res), calculated on cluster (plot) level across 5 subplots
climate_de           <- fread("outData/xy_DE_clim_DWD_1980_2024.csv")
spei_de              <- fread("outData/xy_DE_spei_DWD_1980_2024.csv")

# correcpr precipitation ame
climate <- climate %>% 
  rename(prcp = prec)

# get vegetation data - on subplot level
load("outData/veg.Rdata")


# clean up data -----------------------------------------------------------------
## Vegetation data --------------------------------------
# sum stems first up across vertical groups!!

# create idividual table for each of the vertical regeneration class 
df_stems_saplings <- stem_dens_ha_cluster_sum %>% 
  dplyr::filter(VegType == 'Saplings') %>%  # 20 cm - 2 m
  group_by(cluster) %>% 
  summarise(sum_stems_sapling = sum(total_stems))# %>% 

df_stems_juveniles <- stem_dens_ha_cluster_sum %>% 
  dplyr::filter(VegType == 'Juveniles') %>%  # > 2m heigh
  group_by(cluster) %>% 
  summarise(sum_stems_juvenile = sum(total_stems))# %>% 


df_stems_mature <- stem_dens_ha_cluster_sum %>% 
  dplyr::filter(VegType == 'Mature') %>% 
  group_by(cluster) %>% 
  summarise(sum_stems_mature = sum(total_stems))# %>% 



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

# <<Step 3: Plot - skip, long to run 
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
            drought_prcp = median(prcp, na.rm = T)) 


# merge all predictors together, ID (=subplot) level per year: having a median per years 2018-2023
df_predictors <- 
  climate %>% 
  left_join(climate_subplot_drought, by = join_by(ID)) %>% 
  left_join(dplyr::select(distance_to_edge, c(-country, -dist_ID ))) %>% 
  left_join(dplyr::select(disturbance_chars, c(-region, -country ))) %>% 
  left_join(soil) %>%
  left_join(spei_subplot) %>%          # clim values are medians from 2018-2023
  left_join(spei_subplot_drought) %>%  # clim values are medians from 2018-2020
  left_join(dplyr::select(terrain, c(-country, region, -cluster.x, -cluster.y))) #%>% 
  #mutate(cluster = str_sub(ID, 4, -3)) 



# merge predcitors on cluster level: calculate medians
# keep only temp and prcp: 2021-2023 average
df_predictors_plot <- 
  df_predictors %>% 
  #dplyr::filter(year %in% 2018:2023) %>% 
  group_by(cluster, country) %>% 
  summarise(tmp          = median(tmp, na.rm = T),
            prcp         = median(prcp, na.rm = T),
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



fwrite(df_predictors_plot, 'outData/all_predictors_plot.csv')


# improve climate resolution for germamny and check with original data: ------------
# - get clim anomalies
# - get SPEI
# - calculate medians over 2018-2023

reference_period = 1980:2015

# Calculate temp anomalies
climate_de_anom <- climate_de %>% 
  group_by(cluster) %>%
  mutate(tmp_mean = mean(tmp[year %in% reference_period]),
         tmp_sd = sd(tmp[year %in% reference_period]),
         tmp_z  = (tmp - mean(tmp[year %in% reference_period])) / sd(tmp[year %in% reference_period]),
         prcp_z = (prcp - mean(prcp[year %in% reference_period])) / sd(prcp[year %in% reference_period]))

# get sum of precipitation per year: 
climate_de_prec_sum <- climate_de %>% 
  group_by(cluster, year) %>%
  dplyr::summarize(prcp_sum = sum(prcp, na.rm = T))

# merge precipitation sum back to original table
climate_de_anom_18_23_avg <- climate_de_anom %>% 
  left_join(climate_de_prcp_sum, by = join_by(cluster, year)) %>% 
  dplyr::filter(year %in% 2018:2023) %>% 
  group_by(cluster) %>%
  dplyr::summarise(across(where(is.numeric), \(x) median(x, na.rm = TRUE))) %>%  # Use lambda function for median
  dplyr::select(cluster, tmp, prcp_sum, tmp_z, prcp_z) %>% 
  rename(prcp = prcp_sum)

ggplot(climate_de_anom_18_23_avg, aes(x = prcp, y = tmp)) +
  geom_point() + 
  #geom_jitter() + 
  geom_smooth()
  
# check how many clusters have identical prcp
table(climate_de$prcp)
  

head(spei_de)

spei_de_18_23_avg <- spei_de %>% 
  dplyr::filter(year %in% 2018:2023) %>% 
  ungroup(.) %>% 
  group_by(cluster) %>%
  summarize(spei = median(spei, na.rm  = T))

(spei_de_18_23_avg)

clim_spei_de_upd <- climate_de_anom_18_23_avg %>% 
  full_join(spei_de_18_23_avg) %>% 
  mutate(source = "DWD")


# how often i have clusters with teh same climate values? 
# Count unique clusters where tmp and prcp are identical
identical_clusters_DWD_de <- clim_spei_de_upd %>%
  group_by(tmp, prcp) %>%
  summarise(count = n(), .groups = "drop") %>%  # Count occurrences of each (tmp, prcp) pair
  filter(count > 1)  # Keep only those that appear more than once


identical_clusters_ERA <- df_predictors_plot %>% 
  dplyr::filter(country == '11') %>% 
  group_by(tmp, prcp) %>%
  summarise(count = n(), .groups = "drop") %>%  # Count occurrences of each (tmp, prcp) pair
  filter(count > 1)  # Keep only those that appear more than once

# calculate correlation between tmp and prcp for DWD (grmany) and ERA to make sure that i can replace values;
df_ERA_de <- df_predictors_plot %>% 
  dplyr::filter(country == '11')   %>% 
  mutate(source = 'ERA') %>% 
  dplyr::select(cluster, tmp, prcp) %>% 
  rename(tmp_ERA = tmp,
         prcp_ERA = prcp)

df_de <- clim_spei_de_upd %>% 
  dplyr::select(cluster, tmp, prcp) %>% 
  rename(tmp_DWD = tmp,
         prcp_DWD = prcp) %>% 
  left_join(df_ERA_de) %>% 
  drop_na()
# merge updated data (date with higher climate resolution with original ones): 

# Sample scatter plot
plot(df_de$tmp_DWD, df_de$tmp_ERA, 
     main = "Scatter Plot with Correlation",
     xlab = "tmp_DWD", ylab = "tmp_ERA",
     pch = 19, col = "blue")

# Add a regression line
model <- lm(df_de$tmp_ERA ~ df_de$tmp_DWD)
abline(model, col = "red", lwd = 2)

# Compute correlation
correlation <- cor(df_de$tmp_DWD, df_de$tmp_ERA)

# Add correlation text
text(x = min(df_de$tmp_DWD), y = max(df_de$tmp_ERA), 
     labels = paste("Correlation: ", round(correlation, 2)), 
     pos = 4, col = "red", cex = 1.2)


# Get climate plots for map: TEMP, prcp, SPEI  -------------

### Custom function to return median and IQR for stat_summary
median_iqr <- function(y) {
  median_val <- median(y, na.rm = TRUE)
  iqr_val <- IQR(y, na.rm = TRUE)
  ymin <- median_val - (iqr_val / 2)
  ymax <- median_val + (iqr_val / 2)
  return(c(y = median_val, ymin = ymin, ymax = ymax))
}


# climate full
reference_period <- 1980:2010

# read temp and prcp over years
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
                   prcp = median(prcp)) %>%
  mutate(class = case_when(year %in% 2018:2020 ~ 'drought',
                           TRUE ~ 'ref')) #%>% 

# get reference values
ref_spei <- median(spei_overview$spei[spei_overview$year %in% reference_period], na.rm =T)
ref_tmp  <- median(clim_overview$tmp[clim_overview$year %in% reference_period])
ref_prcp <- median(clim_overview$prcp[clim_overview$year %in% reference_period])


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


# Plot for prcpipitation
p.map.prcp <- clim_overview %>% 
  ggplot(aes(x = year,
             y = prcp,
             color = class)) +
  stat_summary(fun.data = median_iqr, geom = "pointrange", size = 0.1) +  # Use median and IQR for pointrange
  stat_summary(fun = median, geom = "point", size = 0.2) +     # Add point for median
  geom_hline(yintercept = ref_prcp, lty = 'dashed', col = 'grey70') +
  scale_color_manual(values = c('red', 'black')) +
  labs(x = "", y =  expression(paste("Precipitation [mm]", sep=""))) +
  theme_classic() +
  theme(legend.position = "NULL")

windows(7,2)
p.clim.map <- ggarrange(p.map.temp, p.map.prcp,p.map.spei, ncol = 3)
p.clim.map
ggsave(filename = 'outFigs/Fig1_all_clim.png', plot = p.clim.map, width = 7, 
       height = 2, dpi = 300, bg = 'white')



p.clim.map2 <- ggarrange(p.map.temp, p.map.prcp, ncol = 2)
p.clim.map2
ggsave(filename = 'outFigs/Fig1_tmp_prcp.png', plot = p.clim.map2, width = 7, 
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
df_plot_full <- 
  df_plot_veg %>% 
  left_join(df_predictors_plot, by = join_by(cluster)) %>% 
  ungroup() %>% 
  dplyr::select(-c(disturbance_year, 
                   disturbance_agent,
                   silt_extract #, #country,
                   #elevation,
                  # slope, aspect
                  )) %>% 
  dplyr::rename(dominant_species = Species, 
                stem_density     = sum_stems,
                distance_edge    = distance,
                n_vertical       = n_layers) %>% 
  left_join(dat_manag_intensity_cl)# %>% , by = join_by(cluster, management_intensity
  
# set up proper structure (df) fr analysis  and claim characters as factors
df_fin <- as.data.frame(df_plot_full) %>% 
  mutate(cluster = factor(cluster),
         dominant_species = factor(dominant_species),
         country = factor(country))


plot_n <- length(unique(df_plot_full$cluster))

#length(unique(df_plot$cluster))
length(unique(df_plot_full$cluster)) # 849!   - final clusters, 4-5 plots
length(unique(df_predictors_plot$cluster))  # 957 - all clusters, from even with less plots



## add country naming---------------

df_fin <- df_fin %>% 
  rename(prcp = prcp) %>% 
  rename(site = cluster)

# add country indication 
country_regions <- tribble(
  ~country_full, ~region, ~country_abbr,
  "germany", "11, 12, 14, 18, 19, 20, 25", "DE",
  "poland", "17", "PL",
  "czech", "15, 26", "CZ",
  "austria", "13", "AT",
  "slovakia", "16", "SK",
  "slovenia", "23", "SI",
  "italy", "21", "IT",
  "switzerland", "22", "CH",
  "france", "24, 27", "FR"
) %>%
  separate_rows(region, sep = ", ") %>%
  mutate(region = as.integer(region))


df_fin <- df_fin %>%
  mutate(region = as.integer(substr(site, 1, 2)))

# Merge the country information with df_fin
df_fin <- df_fin %>%
  left_join(country_regions, by = "region") %>%
  mutate(country_abbr = case_when(
    site %in% c("24_136", "24_137") ~ "BE",            # Belgium
    site %in% c("24_133", "24_134", "24_135") ~ "LX",  # Luxembourg
    TRUE ~ country_abbr
  )) %>% 
  mutate(country_pooled = case_when( country_abbr == "BE" ~ "FR",  # create pooled data for the eco analysis
                                     country_abbr == "LX" ~ "FR",
                                     TRUE~country_abbr))


fwrite( df_fin, 'outData/indicators_for_cluster_analysis.csv')





# PLOTS: clim cluster : ------------------------
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


# Get summary table 3 env clustersL -------------------------

## ENV CLIM conditions ------------------------------
# Summarize the data based on 'clim_cluster_spei3' for the specified variables

summary_table <- df_fin %>%
  group_by(clim_class) %>%
  summarise(
    temperature = paste0(round(median(tmp, na.rm = TRUE), 2), " (IQR: ", round(IQR(tmp, na.rm = TRUE), 2), ")"),
    precipitation = paste0(round(median(prcp, na.rm = TRUE), 2), " (IQR: ", round(IQR(prcp, na.rm = TRUE), 2), ")"),
    spei = paste0(round(median(spei3, na.rm = TRUE), 2), " (IQR: ", round(IQR(spei3, na.rm = TRUE), 2), ")"),
    #drought_spei = paste0(round(median(drought_spei3, na.rm = TRUE), 2), " (IQR: ", round(IQR(drought_spei3, na.rm = TRUE), 2), ")"),
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


# summary table per drought year (2018):
#summary_table_drought <- 
  df_fin %>%
 # group_by(clim_class) %>%
  summarise(
    drought_temperature = median(drought_tmp, na.rm = TRUE),
    drought_precipitation = median(drought_prcp, na.rm = TRUE),
    drought_spei = median(drought_spei3, na.rm = TRUE)
    ) %>% 
  # ad new columns as reference values 1980-2010
  mutate( ref_tmp = ref_tmp, 
          ref_prcp = ref_prcp,
          ref_spei = ref_spei) %>% 
    mutate(
     temp_diff = drought_temperature - ref_tmp,
      prcp_diff = drought_precipitation - ref_prcp,
      spei_diff = drought_spei - ref_spei
    )

# Print the summary table
print(summary_table_drought)


sjPlot::tab_df(summary_table,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = F,
               file="outTable/clim_cluster_summary_env_conditions.doc",
               digits = 1) 






# 
summary_table_country <- df_fin %>%
group_by(clim_class, country_abbr) %>%
  summarise(
    temperature = paste0(round(median(tmp, na.rm = TRUE), 2), " (IQR: ", round(IQR(tmp, na.rm = TRUE), 2), ")"),
    precipitation = paste0(round(median(prcp, na.rm = TRUE), 2), " (IQR: ", round(IQR(prcp, na.rm = TRUE), 2), ")"),
    spei = paste0(round(median(spei3, na.rm = TRUE), 2), " (IQR: ", round(IQR(spei3, na.rm = TRUE), 2), ")"),
    sand = paste0(round(median(sand_extract, na.rm = TRUE), 2), " (IQR: ", round(IQR(sand_extract, na.rm = TRUE), 2), ")"),
    clay = paste0(round(median(clay_extract, na.rm = TRUE), 2), " (IQR: ", round(IQR(clay_extract, na.rm = TRUE), 2), ")"),
    depth = paste0(round(median(depth_extract, na.rm = TRUE), 2), " (IQR: ", round(IQR(depth_extract, na.rm = TRUE), 2), ")"),
    av_nitrogen = paste0(round(median(av.nitro, na.rm = TRUE), 2), " (IQR: ", round(IQR(av.nitro, na.rm = TRUE), 2), ")")
  ) %>% 
  arrange(country_abbr)

# Print the summary table
print(summary_table_country)

sjPlot::tab_df(summary_table_country,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = F,
               file="outTable/clim_cluster_summary_env_conditions_country.doc",
               digits = 1) 



### Veg indicators Quantiles - no clim cluster -------------------------------------------------
# : summary per quantiles
# Represent results using quantiles, as they are skewed?
qntils = c(0, 0.01, 0.25, 0.5, 0.75, 0.90, 1)

qs_dat_fin <- 
  df_fin %>% 
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
  df_fin %>% 
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
               file="outTable/summary_out_no_clim_cluster.doc",
               digits = 1) 




## Clim: cluster: Veg indicators  -----------------------------------------------

# Table for Mean ± SD
mean_sd_table <- df_fin %>%
  group_by(clim_class) %>%
  summarise(
    stem_density    = paste0(round(mean(stem_density, na.rm = TRUE), 0), " ± ", round(sd(stem_density, na.rm = TRUE), 0)),
    n_vertical      = paste0(round(mean(n_vertical, na.rm = TRUE), 1), " ± ", round(sd(n_vertical, na.rm = TRUE), 1)),
    rIVI            = paste0(round(mean(rIVI, na.rm = TRUE), 1), " ± ", round(sd(rIVI, na.rm = TRUE), 1)),
    richness        = paste0(round(mean(richness, na.rm = TRUE), 0), " ± ", round(sd(richness, na.rm = TRUE), 0))
  )

# Print Mean ± SD table
print(mean_sd_table)


# Table for Median & IQR
median_iqr_table <- df_fin %>%
  group_by(clim_class) %>%
  summarise(
    stem_density    = paste0(round(median(stem_density, na.rm = TRUE), 0), " (IQR: ", round(IQR(stem_density, na.rm = TRUE), 0), ")"),
    n_vertical      = paste0(round(median(n_vertical, na.rm = TRUE), 1), " (IQR: ", round(IQR(n_vertical, na.rm = TRUE), 1), ")"),
    rIVI            = paste0(round(median(rIVI, na.rm = TRUE), 1), " (IQR: ", round(IQR(rIVI, na.rm = TRUE), 1), ")"),
    richness        = paste0(round(median(richness, na.rm = TRUE), 0), " (IQR: ", round(IQR(richness, na.rm = TRUE), 0), ")")
  )

# Print Median & IQR table
print(median_iqr_table)

# for countries

# Table for Mean ± SD
mean_sd_table_country <- df_fin %>%
  group_by(clim_cluster_spei3, country_abbr) %>%
  summarise(
    stem_density    = paste0(round(mean(stem_density, na.rm = TRUE), 0), " ± ", round(sd(stem_density, na.rm = TRUE), 0)),
    n_vertical      = paste0(round(mean(n_vertical, na.rm = TRUE), 1), " ± ", round(sd(n_vertical, na.rm = TRUE), 1)),
    rIVI            = paste0(round(mean(rIVI, na.rm = TRUE), 1), " ± ", round(sd(rIVI, na.rm = TRUE), 1)),
    richness        = paste0(round(mean(richness, na.rm = TRUE), 0), " ± ", round(sd(richness, na.rm = TRUE), 0))
  ) %>% 
  arrange(country_abbr)

# Print Mean ± SD table
print(mean_sd_table_country)


# Table for Median & IQR
median_iqr_table_country <- df_fin %>%
  group_by(clim_class, country_abbr) %>%
  summarise(
    stem_density    = paste0(round(median(stem_density, na.rm = TRUE), 0), " (IQR: ", round(IQR(stem_density, na.rm = TRUE), 0), ")"),
    n_vertical      = paste0(round(median(n_vertical, na.rm = TRUE), 1), " (IQR: ", round(IQR(n_vertical, na.rm = TRUE), 1), ")"),
    rIVI            = paste0(round(median(rIVI, na.rm = TRUE), 1), " (IQR: ", round(IQR(rIVI, na.rm = TRUE), 1), ")"),
    richness        = paste0(round(median(richness, na.rm = TRUE), 0), " (IQR: ", round(IQR(richness, na.rm = TRUE), 0), ")")
  ) %>% 
  arrange(country_abbr)


# Print Median & IQR table
print(median_iqr_table_country)

### export summary tables -------
sjPlot::tab_df(mean_sd_table,
              show.rownames = F,
               file="outTable/clim_cluster_mean_sd_table.doc",
               digits = 2) 

sjPlot::tab_df(mean_sd_table_country,
              show.rownames = F,
               file="outTable/clim_cluster_mean_sd_table_country.doc",
               digits = 2) 

sjPlot::tab_df(median_iqr_table,
               show.rownames = F,
               file="outTable/clim_cluster_median_iqr_table.doc",
               digits = 2) 

sjPlot::tab_df(median_iqr_table_country,
               show.rownames = F,
               file="outTable/clim_cluster_median_iqr_table_country.doc",
               digits = 2) 


