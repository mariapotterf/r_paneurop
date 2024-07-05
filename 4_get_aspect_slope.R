# get slope, aspect

# ------------------------------------------------------------------------------
#   Get terrain characteristics
# ------------------------------------------------------------------------------

library(elevatr)
#library(sf)
library(terra)
library(dplyr)
# define paths
elev_path          <- "rawData/dem"



# list all countries
country_names <- list( "austria","belgium", "czechia", "france", "germany", "italy", "luxembourg", "poland", 
                       "slovakia", "slovenia", "switzerland") #austria" 



get_terrain_chars <- function(country_name) {
  
  print(paste("Processing", country_name))
  #country_name = 'slovakia'
  country <- vect(paste0('outData/dat_', country_name, '.gpkg'))
  xy_proj <- terra::project(country, "EPSG:3035")
  #xy_sf <- st_as_sf(country, "EPSG:3035")
  
  # elevation
  elev_name      <- paste0(toupper(substr(country_name, 1, 1)), tolower(substr(country_name, 2, nchar(country_name))))
  elevation      <- rast(paste0(elev_path, '/dem_', elev_name, '.tif'))
  
  # Check if the current projection of 'elevation' is not 'EPSG:3035'
  if (terra::crs(elevation) != "EPSG:3035") {
    # If it's not, then project it to 'EPSG:3035'
    elevation_proj <- terra::project(x = elevation, "EPSG:3035", method="near")
  } else {
    # If it's already in 'EPSG:3035', just use it as is
    elevation_proj <- elevation
  } 
  
  # Calculate aspect, topography, slope, TRI = terrain roughness index
  slope     <- terra::terrain(elevation, 'slope',     neighbors = 8)
  aspect    <- terra::terrain(slope,     'aspect',    neighbors = 8)
 
  # Extract values to XY points, add IDs
  df_slope     <- as.data.frame(terra::extract(slope, xy_proj, bind = TRUE)) # 
  df_aspect    <- as.data.frame(terra::extract(aspect, xy_proj, bind = TRUE))
  
  # merge all vectors as a new column to XY data
  df <- df_slope %>% 
    left_join(df_aspect, by = join_by(ID, region, country)) 
    
  return(df)
  
  
}

test <- get_terrain_chars('france')
View(test)

out.ls <- lapply(country_names, get_terrain_chars)

# Combine results from all countries
final_terrain <- do.call(rbind, out.ls)

View(final_terrain)
final_terrain %>% 
  dplyr::filter(country == 19 & region == 24) %>% 
  nrow()

# there is 5 NAs in Luxemburg, s they are too close to the border (not enought pixels for 8 neighbors).
# fill in values by their median value:
final_terrain <- final_terrain %>%
  mutate(slope = ifelse(is.na(slope), median(slope, na.rm = TRUE), slope),
         aspect = ifelse(is.na(aspect), median(aspect, na.rm = TRUE), aspect))
 
# total rows: 4755 
# Export data -------------------------------------------------------------
fwrite(final_terrain, 'outData/terrain.csv')

