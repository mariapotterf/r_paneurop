
# Share data : Bc Thesis Benedikt Srb
# Germany

# get stem density: juveniles and saplings (only regeneration layer)
# distance to edge
# DE gpkg
# 

library(data.table)
library(sf)
library(terra)
library(dplyr)

# get vegetation data
load("outData/veg.Rdata")

xy <- vect("outData/xy_clim_cluster.gpkg")
plot(xy)

# filter XY only to germany:  ---------------------------------
xy_de <- xy[xy$country == "DE", ]

# Select specific columns
xy_de <- xy_de[, c("disturbance_severity","distance_edge", "site", "country")]

# View the first few rows
head(xy_de)
nrow(xy_de)

hist(xy_de$distance_edge)
crs(xy_de)

# Convert to EPSG:3035
xy_de_3035 <- terra::project(xy_de, "EPSG:3035")

plot(xy_de_3035)

# filter stem density data ------------------------------------------------------
df_de <- stem_dens_species_long_cluster %>% 
  dplyr::filter(country == "DE") %>% 
  #dplyr::filter(stem_density>0) %>%  # i want to keep in empty sites as well! 
  dplyr::filter(VegType != "Mature") %>% 
  dplyr::rename(site = cluster) %>% 
  dplyr::select(site, 
                VegType, 
                Species, stem_density)

# merge df with xy
merged_spatVect <- merge(xy_de_3035, df_de, by = "site", all.x = TRUE)

# Write to GeoPackage
writeVector(merged_spatVect, 
            'data_share/DE_regen_stem_density.gpkg', 
            filetype = "GPKG",
            overwrite = TRUE)


# export as csv as well ------------------

# Convert SpatVector to a DataFrame
merged_df <- as.data.frame(merged_spatVect)

# Define output CSV file path
csv_output_file <- "data_share/DE_regen_stem_density.csv"

# Export DataFrame to CSV
write.csv(merged_df, file = csv_output_file, row.names = FALSE)


# get number of sites that heve no rereneration --------------------------------
# Filter rows where stem_density is 0
sites_with_only_zero_density <- merged_spatVect %>%
  as.data.frame() %>%          # Convert SpatVector to a data frame for dplyr operations
  group_by(site) %>%           # Group by site
  dplyr::filter(all(stem_density == 0)) %>%  # Keep only sites where all stem_density == 0
  distinct(site)  
# Check the result
length(sites_with_only_zero_density$site)

