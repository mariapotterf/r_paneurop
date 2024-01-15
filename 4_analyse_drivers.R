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
#df_cluster <-
  df_IVI %>% 
  left_join(df_richness) %>%
  left_join(df_vert) %>%
  left_join(stem_dens_ha_cluster)# %>%
  
  
            