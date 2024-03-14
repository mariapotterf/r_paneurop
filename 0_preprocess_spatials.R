

# Prepare spatial data

# load all gpkgs
# load shps of countries
# split shps by country names
# export gpkgs by country

library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(data.table)
library(dplyr)



# list all field data collected; get paths
all_gpkg <- list.files(path = "rawData/collected", 
                       pattern = "*.gpkg$", 
                       full.names = TRUE)  # read the full path to read the files properly

(all_gpkg)

# load gpkgs
all_gpkg.ls <- lapply(all_gpkg, function(name) {vect(name)})


#lapply(all_gpkg, function(name) {crs(name)})

crs(all_gpkg.ls[[12]])


# merge them all in one file
merged_gpkg <- do.call("rbind", all_gpkg.ls)

# convert to sf
merged_gpkg <- sf::st_as_sf(merged_gpkg)


# create an ID and change again to spatial format
# keep only IS and country indetifiers
data_clean <- 
  merged_gpkg %>% 
  dplyr::select(country, region, group, point) %>% 
  ungroup(.) %>% 
  # correct country coding (missing for region 18, which is germany (11))
  mutate(country = case_when(
    region %in% c("11", "12", "14", "18", "19", "20", "25") ~ "11",  # Germany
    region == "17" ~ "12",  # Poland
    region %in% c("15", "26") ~ "13",  # Czech
    region == "13" ~ "14",  # Austria
    region == "16" ~ "15",  # Slovakia
    region == "23" ~ "16",  # Slovenia
    region == "21" ~ "17",  # Italy
    region == "22" ~ "18",  # Switzerland
    region %in% c("24", "27") ~ "19",  # France
    TRUE ~ "NO"  # Default case for any other region
  )) %>%
  dplyr::mutate(group = group + 100) %>% 
  dplyr::mutate(ID = paste(country, region, group, point, sep="_")) %>%  
  #dplyr::mutate(ID_reg = paste(country, region, sep="_") %>% #dplyr::select(5:323) %>% names()
  dplyr::select(ID, country, region)


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
lu <- ne_countries(country = "luxembourg",  type = "countries", returnclass = 'sf', scale = 'medium')
be <- ne_countries(country = "belgium",     type = "countries", returnclass = 'sf', scale = 'medium')

# rename czechia:
cz$sovereignt[cz$sovereignt == "Czech Republic"] <- "Czechia"

# list cuntrues for a lapply loop
cntrs_ls <- list(ge, at, sk, sl,fr, ch, pl, cz, it, lu, be) # 






# Loop over each country 
# function: transform crs system,  clip on country, and export the gpkg data
clip_dat <- function(x, ...) {
  
  # change projection
  x_proj <-st_transform(x, crs(data_clean)) 
  
  # get country name
  name = tolower(x$sovereignt) 
  print(name)
  
  # split the gpkgs by countries, export as individual gpkg per country
  data_clipped <- st_intersection(x_proj, data_clean)
  
  data_clipped <- data_clipped %>% 
    dplyr::select(ID, country, region)
  # export as gpkg
  outName = paste0('dat_', name,  '.gpkg')
  st_write(data_clipped, paste('outData', outName, sep = "/"),
           layer = 'dat', append=FALSE)
  
 
}

# test function -----------------------------------------------------------------
germ_test <- clip_dat(cntrs_ls[[1]])

# export gpkgs by individual countries
lapply(cntrs_ls, clip_dat)








