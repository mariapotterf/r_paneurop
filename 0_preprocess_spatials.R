

# Prepare spatial data

# load all gpkgs
# load shps of countries
# split shps by country names
# export gpkgs by country

library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)


setwd('C:/Users/ge45lep/Documents/2023_PanEuropean/r_paneurop/rawData/collected')


# list all field data collected; get paths
all_gpkg <- list.files(pattern = "*.gpkg")

# load gpkgs
all_gpkg.ls <- lapply(all_gpkg, function(name) {vect(name)})

# merge them all in one file
merged_gpkg <- do.call("rbind", all_gpkg.ls)

# set as sf
merged_gpkg <- sf::st_as_sf(merged_gpkg)


# get shp of european countries
ge <- ne_countries(country = "germany",     type = "countries", returnclass = 'sf', scale = 'medium')
at <- ne_countries(country = "austria",     type = "countries", returnclass = 'sf', scale = 'medium')
sk <- ne_countries(country = "slovakia",    type = "countries", returnclass = 'sf', scale = 'medium')
sl <- ne_countries(country = "slovenia",    type = "countries", returnclass = 'sf', scale = 'medium')
fr <- ne_countries(country = "france",      type = "countries", returnclass = 'sf', scale = 'medium')
ch <- ne_countries(country = "switzerland", type = "countries", returnclass = 'sf', scale = 'medium')
pl <- ne_countries(country = "poland",      type = "countries", returnclass = 'sf', scale = 'medium')
cz <- ne_countries(country = "czech republic", type = "countries", returnclass = 'sf', scale = 'medium')
it <- ne_countries(country = "italy",       type = "countries", returnclass = 'sf', scale = 'medium')

# merge EU countries
mrg_cntrs <- dplyr::bind_rows(list(ge, at, sk, sl,fr, ch, pl, cz, it))

# project vector data to the same crs
mrg_cntrs_pr <-st_transform(mrg_cntrs, crs(merged_gpkg)) 


plot(mrg_cntrs_pr$geometry)
plot(merged_gpkg$geometry, add = T)


# FITS nicely!
  # 

# split the gpkgs by countries, export




