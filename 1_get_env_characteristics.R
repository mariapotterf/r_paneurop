
# 
# Prepare spatial database of the each cluster point

# read field data by level of country
# get raster data:
# elevation
# climat vars - ERA-NEt - data already downloaded
# data from the Cornelius: disturbances
#    - year
#    - agent
#    - distance to nearest edge (merge 2018-2020 into one!)


library(terra)

# list all countries
country_name = c('austria')

# list all countries
countries_ls <- list( "austria", "czechia", "france", "germany", "italy", "poland", 
                "slovakia", "slovenia", "switzerland") #austria" 


# read the plots by country
# read respective rasters: DEM, disturbances...
#    = stack them
# extract info by plot level
# move toanother country

# get paths to rasters -----------------------------------------
elev_path          <- "rawData/dem"
#dist_agent_path    <- "rawData/disturb_data"
dist_path          <- "rawData/disturb_data" # for year and severity


# read field data for one country --------------------------------
country = vect(paste0('outData/dat_', country_name, '.gpkg'))

# read raster data

# disturbance year
disturb_name = paste0('disturbance_year_', country_name, '.tif')
disturbance  = rast(paste(dist_path, country_name, disturb_name, sep = '/'))

# disturbace severity
severity_name = paste0('disturbance_severity_', country_name, '.tif')
severity      = rast(paste(dist_path, country_name, severity_name, sep = '/'))

# add agent !!!
# elevation
elev_name      <- paste0(toupper(substr(country_name, 1, 1)), tolower(substr(country_name, 2, nchar(country_name))))
elevation      <- rast(paste0(elev_path, '/dem_', elev_name, '.tif'))
elevation_proj <- project(x = elevation, y = disturbance, method="near")
elevation_proj <- resample(x = elevation_proj, y = disturbance, method="near")


desired_crs <- crs(disturbance)

if (crs(severity) != desired_crs) {
  severity <- terra::project(severity, desired_crs)
}

if (crs(elevation_proj) != desired_crs) {
  elevation_proj <- terra::project(elevation_proj, desired_crs)
}








# Process data ------------------------------------------------------
country_proj <- project(country, "EPSG:3035")

# create raster stacks
dist.stack <- c(disturbance, severity, elevation_proj)
names(dist.stack) <- c("disturbance_year", "disturbance_severity", "elevation")


# extract elevation for every point
plots.disturbance <- extract(dist.stack, country_proj, method = "simple", bind=TRUE)






