
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







