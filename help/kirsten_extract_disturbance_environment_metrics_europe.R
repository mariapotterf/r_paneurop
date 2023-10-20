#---------- libaries

library(terra)
library(dplyr)


setwd("path_to_your_plot_points") # monitoring points


# --- load centroids with patch size

plots <- read.csv("plots.csv", header=T, sep=",")
plots <- vect(plots, geom=c("longitude", "latitude"), crs = "EPSG:4326") #create point layer from csv # set correct coordiante names and CRS

# reproject to CRS used with disturbance map: 
#// depending in which CRS your plots and the other environmental data are, you can also think about reprojecting the disturbance maps

plots <- project(plots, "EPSG:3035")


# set paths to disturbance maps and country shapes

path_dist_maps <- "F:/Projects/DisturbanceMappingEurope/mapping/results/version1.1/" # to be clean, I would rather copy-paste this folder to your homefolder (maybe ask Cornelius first)

path_countries <- "F:/Projects/DisturbanceMappingEurope/mapping/data/gis/countries/" # I would copy paste this folder to your local home, super useful to have al the country shapes!


#countries <- list.files(path_dist_maps) #if you want to create a vector with all country names

# subset to countries with observation plots:
countries <- c(  "belgium", "czechia", "france", "germany", "hungary", "ireland", "netherlands", "norway", "poland", "slovakia", "slovenia", 
                 "sweden","switzerland","ukraine","unitedkingdom") #austria" 


# initialize the first extraction as input for the loop

cntr <- "austria"

disturbance <- rast(paste0(path_dist_maps, cntr, "/disturbance_year_", cntr, ".tif"))
severity <- rast(paste0(path_dist_maps, cntr, "/disturbance_severity_", cntr, ".tif"))

dist.stack <- c(disturbance, severity, agent)
names(dist.stack) <- c("disturbance_year", "disturbance_severity", "disturbance_agent")

border <- vect(paste0(path_countries, cntr, ".shp"))
plots.sub <- mask(plots, border)

plots.disturbance <- extract(dist.stack, plots.sub, method = "simple", bind=TRUE)


# --- extract disturbance infos for the other countries


for (i in 1:length(countries)) {
  
  cntr <- countries[i]
  
  print(cntr)
  
  disturbance <- rast(paste0(path_dist_maps, cntr, "/disturbance_year_", cntr, ".tif"))
  severity <- rast(paste0(path_dist_maps, cntr, "/disturbance_severity_", cntr, ".tif"))
  # load all the data to extract
  
  border <- vect(paste0(path_countries, cntr, ".shp")) # load border of the country
  
  #mask all data layers to the country
  disturbance <- mask(disturbance, border) 
  #.....
  
  #stack the data 
  dist.stack <- c(disturbance, severity)
  names(dist.stack) <- c("disturbance_year", "disturbance_severity")  # add variabel names how you want them in the final csv
  
  plots.sub <- mask(plots, border) # mask your european-wide plots to the specific country
  
  
  plots.sub <- extract(dist.stack, plots.sub, method = "simple", bind=TRUE) # extract raster values for each plot point coordinate and bind them as columns to plot df
  
  plots.disturbance <- rbind(plots.disturbance, plots.sub) # bind your sub df for the current country to oyur previous created df for austria (with same data structure!)
}
#}

df.dist.points <- as.data.frame(plots.disturbance) # create a data.frame out of the points

# store your data frame with the information