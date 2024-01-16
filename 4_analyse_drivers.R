# create final database
# merge veg data
# predictors - climate, disturbance chars, soils
# dependent - regeneretion, density, richness
# output file
# 

library(data.table)
library(dplyr)
library(terra)
library(ggplot2)

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

# cluster level
df_cluster <-
  df_IVI %>% 
  left_join(df_richness) %>%
  left_join(df_vert) %>%
  left_join(stem_dens_ha_cluster)# %>%


# Analysis ------------------------------------------------------------
# 1. How many clusters are only manag/unmanag/mixed?
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

#Managed    Mixed Unmanaged 
# 668       155        84 

## 2. species composition per stems: saplings vs juveniles?? --------------
# how many plots with Mature trees?
# 
dom_species <- stem_dens_species_long %>% 
  filter(VegType != 'Survivor' ) %>%  # Survivors: on 191
  group_by(Species, VegType, manag) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup(.) %>% 
  group_by(VegType) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share = sum_stems/sum_vegType*100) #%>% 
  
dom_species %>% 
  dplyr::filter(share > 1) %>% 
  ggplot(aes(x = VegType,
             y = share,
             fill = Species)) +
  geom_col()
  
# see all species
dom_species %>% 
  #dplyr::filter(share > 1) %>% 
  ggplot(aes(x = Species,
             y = share,
             fill = VegType)) +
  geom_bar(position = "fill", stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #geom_col('')



# all species by aboundance
dom_species %>% 
  #dplyr::filter(share > 1) %>% 
  ggplot(aes(x = Species,
             y = share,
             fill = VegType)) +
  geom_bar(stat = "identity")+
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


# 3. how many clusters/plots has less then 2000 regen/ha?
# which country is the most affected? (where they ccur most often?)


# 4. stems vs weather - mean 2019-2022, or SPEI

# 5. management effect: compare manag vs unmanaged stem densities
  
            