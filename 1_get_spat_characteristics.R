
# 
# Prepare spatial database of the each cluster point

# get data:
# elevation
# climat vars - ERA-NEt - data already downloaded
# data from the Cornelius: disturbances
#    - year
#    - agent
#    - distance to nearest edge (merge 2028-2020 into one!)


# how to run this effectively on server?
# access to rasters 
#
# how to share my script?
# where is my output data?
# paths issues?

# to do:
# doeanload all data
# split the gpkg points by country
# loop everything on the level of country or th eminimal bounding box
# get shps from Cornelius



library(terra)

setwd('C:/Users/ge45lep/Documents/2023_PanEuropean/r_paneurop')
#setwd('/outData')

all_gpkg <- list.files(path = "outData", pattern = "*.gpkg$")
# list all field data collected; get paths
all_gpkg <- list.files(pattern = "*.gpkg$")
(all_gpkg)
# load gpkgs, if there is no empty geometry
#all_gpkg.ls <- lapply(all_gpkg, function(name) {vect(name)})

all_gpkg.ls <- lapply(all_gpkg, function(name) {
  # Try to read the file
  tryCatch({
    vect_data <- vect(name)
    # Return the data if successful
    return(vect_data)
  }, error = function(e) {
    # If there's an error, return NULL or another appropriate value
    return(NULL)
  })
})

# Filter out NULL entries or whatever value you decided to return on error
all_gpkg.ls <- Filter(Negate(is.null), all_gpkg.ls)

# project collected field data
plots_ls_proj <- lapply(all_gpkg.ls, function(name) {project(name, "EPSG:3035")})
#plots <- project(plots, "EPSG:3035")

# get elevation, disturbance severity and agent


# read elevation maps --------------------------------------------------

ras_path <- "rawData/dem/dem_Austria.tif"
ras <- terra::rast(ras_path)

#countries <- list.files(path_dist_maps) #if you want to create a vector with all country names

# subset to countries with observation plots:
countries <- c( "Austria", "Czechia", "France", "Germany", "Italy", "Poland", 
                "Slovakia", "Slovenia", "Switzerland") #austria" 


# initialize the first extraction as input for the loop

cntr <- "austria"

disturbance <- rast(paste0(path_dist_maps, cntr, "/disturbance_year_", cntr, ".tif"))
severity <- rast(paste0(path_dist_maps, cntr, "/disturbance_severity_", cntr, ".tif"))

dist.stack <- c(disturbance, severity, agent)
names(dist.stack) <- c("disturbance_year", "disturbance_severity", "disturbance_agent")

border <- vect(paste0(path_countries, cntr, ".shp"))
plots.sub <- mask(plots, border)

plots.disturbance <- extract(dist.stack, plots.sub, method = "simple", bind=TRUE)




