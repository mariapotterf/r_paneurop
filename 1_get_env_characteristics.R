
# 
# Prepare spatial database of the each cluster point
# 

# read field data by level of country
# get raster data:
# elevation
# climat vars - ERA-NEt - data already downloaded
# data from the Cornelius: disturbances
#    - year
#    - agent - NA = no disturbance
# 1 = bark beetle or wind disturbances (both classes had to be grouped due to technical reasons)
# 2 = fire disturbances
# 3 = other disturbances, mostly harvest but might include salvage logging go small-scale natural disturbances and infrequent other natural agents (e.g., defoliation, avalanches, etc.)


#    - distance to nearest edge ? (merge 2018-2020 into one!)


library(terra)
library(raster)


# list all countries
country_name = c('austria')















# read the plots by country
# read respective rasters: DEM, disturbances...
#    = stack them
# extract info by plot level
# move toanother country

# get paths to rasters -----------------------------------------
elev_path          <- "rawData/dem"
dist_path          <- "rawData/disturb_data" # for year and severity


# function to loop over ----------------------------------------


# list all countries
country_names <- list( "austria", "czechia", "france", "germany", "italy", "poland", 
                       "slovakia", "slovenia", "switzerland") #austria" 






# distances will contain the distance of each point to the nearest edge of its patch
library(sf)

library(terra)
library(raster)



# Create a small example raster
mat <- matrix(c(1, NA, 2, 2, 3,
                NA, 2, 2, 2, 4,
                2, 2, 2, 2, 2,
                2, 2, 2, 1, 2,
                2, 2, 2, 2, 2), byrow = T, nrow=5, ncol=5)
example_raster <- rast(nrows=5, ncols=5, vals=mat, crs = 'EPSG:3035')

# Reclassify the raster - change values '2' to NA
reclass_matrix <- matrix(c(1,2, NA), ncol=3, byrow=TRUE)  #matrix(c(1, NA), ncol=2, byrow=TRUE)
reclassified_raster <- classify(example_raster, reclass_matrix)

# Assuming reclassified_raster is your raster
# Create a data frame with the coordinates of the point
point_df <- data.frame(x = 150, y =-80)

# Create a SpatVector from the data frame
point_vector <- vect(point_df, geom = c("x", "y"), crs = crs(reclassified_raster))


# Calculate the distance to the nearest NA
distance_to_na <- distance(reclassified_raster) # , point_vector


# Plotting to visualize the results
plot(example_raster, main="example_raster")
plot(reclassified_raster, main="Reclassified Raster")
plot(distance_to_na, main="")
plot(point_vector, main="", add = T)

# Measure the distance from the point to the nearest NA
point_distance <- extract(distance_to_na, point_vector)
(point_distance)
plot(distance_to_na, main="Distance to Nearest NA")
plot(point_vector, main="", add = T)









#  example with lines

library(terra)
set.seed(123)

# Create an empty raster and set all values to NA initially
raster <- rast(ext=c(401000, 402000, 4405000, 4406000), res = 30, crs="local")
values(raster) <- NA

# Define a patch area and set these cells to 2 (covering about 1/4 of the raster)
patch_cells <- cells(raster, ext(c(401000, 401500, 4405000, 4405500)))
raster[patch_cells] <- 2

# Create points that fall within the patch
points <- data.frame(
  ID = 1:4,
  x = runif(4, 401000, 401500),  # X coordinates within the patch
  y = runif(4, 4405000, 4405500) # Y coordinates within the patch
)

points <- data.frame(
  ID = 1,
  x = 401200,  # X coordinates within the patch
  y = 4405100 # Y coordinates within the patch
)

pts <- vect(points, geom=c("x", "y"), crs=crs(raster))

# Plotting to visualize the raster and points
plot(raster, main="Raster with Patch")
plot(lns, add = TRUE, col = 'blue', lwd = 4)
plot(pts, add=TRUE, col="red", pch=16)


# convert raster patches to polygons and to lines to measure the distances
lns <- as.lines(as.polygons(raster, dissolve=TRUE))

nearest(pts, lns) |> values()

x <- distance(raster>0)
extract(x, pts)













# distance to edge: test for slovakia -----------------------------------------------------------
country_name = 'slovakia'

print(country_name)

# read field data 
country = vect(paste0('outData/dat_', country_name, '.gpkg'))
country_proj <- project(country, "EPSG:3035")

# disturbance year
disturb_name = paste0('disturbance_year_', country_name, '.tif')
disturbance  = rast(paste(dist_path, country_name, disturb_name, sep = '/'))

# read raster data
desired_crs <- crs(disturbance)

# Create a SpatVector from the data frame
point_vector <- country_proj[1,]# vect(point_df, geom = c("x", "y"), crs = crs(reclassified_raster))


# get buffer
buff <- buffer(point_vector, 100)

# crop data and then mask them to have a circle
disturbance_crop <- crop(disturbance, buff)
disturbance_mask <- mask(disturbance_crop, buff)

plot(disturbance_mask)

# Reclassify the raster - change values '2' to NA
reclass_matrix <- matrix(c(1985, 2017.1, 1,
                           2017.5, 2021, NA), ncol=3, byrow=TRUE)
reclassified_raster <- classify(disturbance_mask, reclass_matrix)

plot(reclassified_raster)

# Calculate the distance to the nearest NA
distance_to_na <- distance(reclassified_raster)

# Plotting to visualize the results
plot(reclassified_raster, main="Reclassified Raster")
plot(distance_to_na, main="Distance to Nearest NA")

# Measure the distance from the point to the nearest NA
point_distance <- extract(distance_to_na, point_vector)
(point_distance)
plot(distance_to_na, main="Distance to Nearest NA")
plot(point_vector, main="", add = T)





# Make a function distance ti edge -------------------------
process_point <- function(point, disturbance, buffer_dist) {
  # Create buffer around the point
  buff <- buffer(point, buffer_dist)
  
  # Crop and mask the disturbance raster
  disturbance_crop <- crop(disturbance, buff)
  disturbance_mask <- mask(disturbance_crop, buff)
  
  # Reclassify the raster
  reclass_matrix <- matrix(c(1985, 2017.1, 1,
                             2017.5, 2021, NA), ncol=3, byrow=TRUE)
  reclassified_raster <- classify(disturbance_mask, reclass_matrix)
  
  # Calculate the distance to the nearest NA
  distance_to_na <- distance(reclassified_raster)
  
  # Extract the distance for the point
  point_distance <- extract(distance_to_na, point)
  
  return(data.frame(ID = point$ID, distance = point_distance))
}

# Example usage:
# Assuming 'disturbance' is your raster and 'country_proj' is your SpatVector

results <- lapply(1:nrow(country_proj), function(i) {
  point_vector <- country_proj[i,]
  process_point(point_vector, disturbance, 500)
})

# Combine all results into one dataframe
final_results <- do.call(rbind, results)
print(final_results)




# Assuming reclassified_raster is your raster
# Create a data frame with the coordinates of the point
point_df <- data.frame(x = 150, y =-80)

# Create a SpatVector from the data frame
point_vector <- vect(point_df, geom = c("x", "y"), crs = crs(reclassified_raster))


# Calculate the distance to the nearest NA
distance_to_na <- distance(reclassified_raster) # , point_vector


# Plotting to visualize the results
plot(example_raster, main="example_raster")
plot(reclassified_raster, main="Reclassified Raster")
plot(distance_to_na, main="")
plot(point_vector, main="", add = T)

# Measure the distance from the point to the nearest NA
point_distance <- extract(distance_to_na, point_vector)
(point_distance)
plot(distance_to_na, main="Distance to Nearest NA")
plot(point_vector, main="", add = T)












# function to extract all data
extract_env_info <- function(country_name) {
   print(country_name)
 
   # read field data 
   country = vect(paste0('outData/dat_', country_name, '.gpkg'))
   country_proj <- project(country, "EPSG:3035")
   
   # read raster data
   desired_crs <- crs(disturbance)
   
   # disturbance year
   disturb_name = paste0('disturbance_year_', country_name, '.tif')
   disturbance  = rast(paste(dist_path, country_name, disturb_name, sep = '/'))
   
   # disturbace severity
   severity_name = paste0('disturbance_severity_', country_name, '.tif')
   severity      = rast(paste(dist_path, country_name, severity_name, sep = '/'))
   
   # disturbance agent
   agent_name = paste0('fire_wind_barkbeetle_', country_name, '.tif')
   agent      = rast(paste(dist_path, agent_name, sep = '/'))
   crs(agent) <-desired_crs
   
   # elevation
   elev_name      <- paste0(toupper(substr(country_name, 1, 1)), tolower(substr(country_name, 2, nchar(country_name))))
   elevation      <- rast(paste0(elev_path, '/dem_', elev_name, '.tif'))
   elevation_proj <- terra::project(x = elevation, y = disturbance,  method="near")
   elevation_proj <- terra::resample(x = elevation_proj, y = disturbance, method="near")
   
   # create raster stacks
   dist.stack <- c(disturbance, severity, agent, elevation_proj)
   names(dist.stack) <- c("disturbance_year", 
                          "disturbance_severity", 
                          "disturbance_agent", 
                          "elevation")
 
    # extract elevation for every point
   plots.disturbance <- extract(dist.stack, country_proj, method = "simple", bind=TRUE)
   
   # export plot
   return(plots.disturbance)
   
 }

out <- extract_env_info('austria')

# list all countries
country_names <- list( "austria", "czechia") #austria" 

out_ls <- lapply(country_names, extract_env_info)

# merge them all in one file
merged_ls <- do.call("rbind", out_ls)

# export final table
save(merged_ls, file="outData/plots_env.Rdata")









# Test for single country -------------------------------------------------


# read field data for one country --------------------------------
country = vect(paste0('outData/dat_', country_name, '.gpkg'))

# read raster data
# make sure the resolution is the same
desired_crs <- crs(disturbance)

# disturbance year
disturb_name = paste0('disturbance_year_', country_name, '.tif')
disturbance  = rast(paste(dist_path, country_name, disturb_name, sep = '/'))

# disturbace severity
severity_name = paste0('disturbance_severity_', country_name, '.tif')
severity      = rast(paste(dist_path, country_name, severity_name, sep = '/'))

agent_name = paste0('fire_wind_barkbeetle_', country_name, '.tif')
agent      = rast(paste(dist_path, agent_name, sep = '/'))
crs(agent) <-desired_crs

elev_name      <- paste0(toupper(substr(country_name, 1, 1)), tolower(substr(country_name, 2, nchar(country_name))))
elevation      <- rast(paste0(elev_path, '/dem_', elev_name, '.tif'))
elevation_proj <- terra::project(x = elevation, y = disturbance,  method="near")
#elevation_proj <- resample(x = elevation_proj, y = disturbance, method="near")



# 
# if (crs(severity) != desired_crs) {
#   severity <- terra::project(severity, desired_crs)
# }
# 
# if (crs(elevation_proj) != desired_crs) {
#   elevation_proj <- terra::project(elevation_proj, desired_crs)
# }

crs(elevation_proj) == crs(disturbance)
crs(elevation_proj) == crs(agent)


# Process data ------------------------------------------------------
country_proj <- project(country, "EPSG:3035")

# create raster stacks
dim(disturbance)
dim(severity)
dim(agent)
dim(disturbance)

crs(disturbance)
crs(severity)
crs(agent)
crs(elevation_proj)

dist.stack <- c(disturbance, severity, agent, elevation_proj)
names(dist.stack) <- c("disturbance_year", 
                       "disturbance_severity", 
                       "disturbance_agent", 
                       "elevation")


# extract elevation for every point
plots.disturbance <- extract(dist.stack, country_proj, method = "simple", bind=TRUE)







