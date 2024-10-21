# Extract future species distribution 
# wessely 2024

# process;

# read the list of all Wessely species (67)
# read my species list (~37)
# narrows down speies list
# read the vect points

# look over the rasters to extract the values
# this will be 0/1 per species and decade
# merge 

library("terra")
library(tidyr)
library(ggplot2)
library(viridis)
library(raster)
library(reshape2)



# Get species Wessely and field data - prepared and merged manually
species<- read.csv("rawData/tree_sp_field_wessely_merged.csv", sep = ';')[,3]

# test single one 
# read locations:
xy <- vect("outData/xy_clim_cluster.gpkg")
xy_proj <- project(xy, crs(r))

# keep only site  colums
xy <- xy_proj[, "site"]  

# plot overlap
plot(r)
plot(xy, add = T)

# get raster naming 
species = c("Abiesalba", "Acercampestre")
## define climate change scenarios
scenarios <- c("rcp26","rcp45","rcp85")
## define timesteps
timesteps <- seq(2020,2090,by=10)
## to delet?

# read file
rastPath <- "rawData/wessely_species_bottleneck/data_submit/data_submit/modelResults"
r <- rast(paste(rastPath, raster_name_in, sep = '/'))


raster_name_in = paste(one_species, '_rcp26_', '2020.tif', sep = '')
rast_name_out  = paste(one_species, '_rcp26_', '2020', sep = '')
# extrcat the raster value to each of teh sites
# Extract raster values at the locations of xy_projected
df <- extract(r, xy, ID = F, 
                            bind = T    # add sites names
)

names(df)[names(df) == "layer"] <- rast_name_out
              
              
              # View the extracted values
head(df)



# test loop


process_raster <- function(species, scenario, timestep, xy) {
  # Create raster input and output names
  raster_name_in <- paste0(species, '_', scenario, '_', timestep, '.tif')
  rast_name_out  <- paste0(species, '_', scenario, '_', timestep)
  
  # Define the raster path
  rastPath <- "rawData/wessely_species_bottleneck/data_submit/data_submit/modelResults"
  
  # Read the raster
  r <- rast(file.path(rastPath, raster_name_in))
  
  # Extract raster values to XY points
  df <- extract(r, xy, ID = FALSE, bind = TRUE)  # 'ID = FALSE' avoids including the raster ID
  
  # Rename the 'layer' column to the raster name
  names(df)[names(df) == "layer"] <- rast_name_out
  
  df <- as.data.frame(df)
  
  # Return the data frame
  return(df)
}

# Define species, scenarios, and timesteps
#species <- c("Abiesalba", "Acercampestre")
scenarios <- c("rcp26", "rcp45", "rcp85")
timesteps <- seq(2020, 2090, by = 10)

# List to hold the extracted data frames
all_data <- list()

# Use nested lapply to loop over all combinations
all_data <- lapply(species, function(sp) {
  lapply(scenarios, function(scenario) {
    lapply(timesteps, function(timestep) {
      # Apply the function to each combination
      process_raster(sp, scenario, timestep, xy)
    })
  })
})

# Combine the list of data frames into one data frame by "site"
combined_df <- Reduce(function(x, y) merge(x, y, by = "site", all = TRUE), all_data_flat)

# View the combined data frame
head(combined_df)

