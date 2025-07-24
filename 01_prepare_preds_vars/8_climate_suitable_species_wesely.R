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

gc()

library("terra")
library(tidyr)
library(ggplot2)
library(viridis)
library(raster)
library(reshape2)
library(data.table)



# Get species Wessely and field data - prepared and merged manually
# 3rd columns has Wessely species
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

# TEST - get raster naming --------------------------------------------------- 
#species = c("Abiesalba", "Acercampestre")
## define climate change scenarios
scenarios <- c("rcp26","rcp45","rcp85")
## define timesteps
timesteps <- seq(2020,2090,by=10)
## to delet?

# read file
rastPath <- "rawData/wessely_species_bottleneck/data_submit/data_submit/modelResults"
r <- rast(paste(rastPath, raster_name_in, sep = '/'))
my_crs <- crs(r)

raster_name_in = paste(one_species, '_rcp26_', '2020.tif', sep = '')
rast_name_out  = paste(one_species, '_rcp26_', '2020', sep = '')
# extrcat the raster value to each of teh sites
# Extract raster values at the locations of xy_projected
df <- extract(r, xy, ID = F, 
                            bind = T    # add sites names
)
# test loop


process_raster <- function(species, scenario, timestep, xy) {
  # Create raster input and output names
  raster_name_in <- paste0(species, '_', scenario, '_', timestep, '.tif')
  rast_name_out  <- paste0(species, '_', scenario, '_', timestep)
  
  # Define the raster path
  rastPath <- "rawData/wessely_species_bottleneck/data_submit/data_submit/modelResults"
  
  # Read the raster
  r <- rast(file.path(rastPath, raster_name_in))
  
  # set crs
  crs(r) <- my_crs
  
  # Extract raster values to XY points
  df <- extract(r, xy, ID = F, bind = TRUE)  # 'ID = FALSE' avoids including the raster ID
  
  # Rename the 'layer' column to the raster name
  names(df)[names(df) == "layer"] <- rast_name_out
  
 # print(rast_name_out)
  # Check if CRS between raster and vector are different
 # if (crs(r) != crs(xy)) {
  #  print("CRS mismatch between raster and vector!")
  #}
  
  df <- as.data.frame(df)
  
  # Return the data frame
  return(df)
}


# Define species, scenarios, and timesteps
species2 <- c("Tiliacordata", "Quercusfrainetto", "Ostryacarpinifolia")
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

# Flatten the nested list structure into a single list of data frames
flattened_data <- do.call(c, do.call(c, all_data))

# Combine all data frames by "site" column
combined_df <- Reduce(function(x, y) merge(x, y, by = "site", all = TRUE), flattened_data)


# View the combined data frame
dim(combined_df)
names(combined_df)


# Convert the wide dataframe to long format
combined_long_df <- combined_df %>%
  pivot_longer(
    cols = -site,                # Keep 'site' as the ID column
    names_to = "raster_name",     # Create a new column for the raster names
    values_to = "raster_value"    # Create a new column for the raster values
  )

# Split the 'raster_name' column into 'species', 'scenario', and 'timestep' using '_'
future_species_full <- combined_long_df %>%
  separate(raster_name, into = c("species", "scenario", "timestep"), sep = "_")

# calculate if teh species is continuous: sum per all years == 8 
# calculate the presence values per site, scenario, species over all years
future_species_lifecycle <- 
  future_species_full %>% 
  group_by(site, species, scenario) %>% 
    summarise(occurence_years = sum(raster_value), .groups = 'drop') %>%
    mutate(overall_presence = case_when(
      occurence_years == 8 ~ 1,
      TRUE ~ 0
    ))

# export final table
fwrite(future_species_lifecycle, 'outTable/species_presence_clim_change.csv')
fwrite(future_species_full, 'outTable/species_presence_clim_change_full.csv')

