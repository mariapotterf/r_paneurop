# create final database
# merge veg data
# predictors - climate, disturbance chars, soils
# dependent - regeneretion, density, richness
# output file
# 

library(data.table)
library(dplyr)
library(terra)

climate           <- fread("outData/climate_18_23.csv")
distance_to_edge  <- fread("outData/distance_to_edge.csv")
disturbance_chars <- fread("outData/disturbance_chars.csv")
soil              <- as.data.frame(terra::vect("rawData/extracted_soil_data/extracted_soil_data.gpkg"))
terrain           <- fread("outData/terrain.csv")


load("outData/veg.Rdata")

disturbance_chars <- disturbance_chars %>%
  mutate(country = as.character(country))

# merge all predictors together, ID level
df <- 
  climate %>% 
  left_join(select(distance_to_edge, c(-country, -dist_ID ))) %>% 
  left_join(select(disturbance_chars, c(-region, -country ))) %>% 
  left_join(soil) %>%
  left_join(select(terrain, c(-country))) %>% 
  left_join(select(stem_dens_species_long, c(-country)))


# select the dominant species per cluster
df_IVI_max <- df_IVI %>% 
  dplyr::select(cluster, Species, manag, country,rIVI) %>% 
  ungroup(.) %>% 
  group_by(cluster, country, manag) %>% 
  filter(rIVI == max(rIVI)) %>% 
  slice(1)

# how many clusters I haver only managed, only unmanaged and how many mixed??
table(df_IVI_max$cluster, df_IVI_max$manag) %>% 
  View()

# cluster level
df_cluster <-
  df_IVI %>% 
  left_join(df_richness) %>%
  left_join(df_vert) %>%
  left_join(stem_dens_ha_cluster)# %>%


# Analysis ------------------------------------------------------------
# 1. How many clusters are only manag/unmanag/mixed?
n_clust_manag <- stem_dens_species_long %>%
  ungroup(.) %>% 
  dplyr::select(cluster, manag) %>% 
  group_by(cluster, manag) %>% 
  distinct(.) %>%
  ungroup(.) %>% 
  group_by(cluster) %>% 
  summarise(n_manag = n()) #%>% 
  #View()
  #filter(n_manag>1)
  
table(n_clust_manag$n_manag)

manag_types <- stem_dens_species_long %>%
  ungroup(.) %>% 
  dplyr::select(cluster, manag) %>% 
  group_by(cluster, manag) %>% 
  distinct(.) %>%
  ungroup(.) %>% 
  group_by(cluster) %>% 
  summarise(type = case_when(
    all(manag == 'Unmanaged') ~ 'Unmanaged',
    all(manag == 'Managed') ~ 'Managed',
    TRUE ~ 'Mixed'
  )) %>%
  ungroup(.) %>% 
  dplyr::select(-cluster) %>% 
  table()

# 907 clusters in total
# single manag type: 752 clusters
# mixed manag type: 155 clusters


## 2. species composition per stems: saplings vs juveniles?? --------------
# how many plots with Mature trees?
# 
stem_dens_species_long %>% 
  filter(VegType != 'Survivor' )  # 191


# 3. how many clusters/plots has less then 2000 regen/ha?
# 4. stems vs weather - mean 2019-2022
  
            