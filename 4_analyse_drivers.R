# create final database
# merge veg data
# predictors - climate, disturbance chars, soils
# dependent - regeneretion, density, richness
# output file
# 

library(data.table)
library(dplyr)
library(terra)

climate           <- fread("outData/climate.csv")
distance_to_edge  <- fread("outData/distance_to_edge.csv")
disturbance_chars <- fread("outData/disturbance_chars.csv")
soil              <- as.data.frame(terra::vect("rawData/extracted_soil_data/extracted_soil_data.gpkg"))
terrain           <- fread("outData/terrain.csv")


disturbance_chars <- disturbance_chars %>%
  mutate(country = as.character(country))

# merge all predictors together
df <- 
  climate %>% 
  left_join(select(distance_to_edge, c(-country, -dist_ID ))) %>% 
  left_join(select(disturbance_chars, c(-region, -country ))) %>% 
  left_join(soil) %>%
  left_join(select(terrain, c(-country)))
  
            