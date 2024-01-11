# get climatic characteristics
# temp
# prec
# temp anomalies???
# from ERA-NET data


library(terra)
library(sf)
library(tidyr)
library(dplyr)
library(data.table)
library(PCICt)
library(zoo)                # for as.Date() specification


# get data from Pan-Europe ----------------------------------------------------


#list all field data collected; get paths
all_gpkg <- list.files(path = "rawData/collected", 
                       pattern = "*.gpkg$", 
                       full.names = TRUE)  # read the full path to read the files properly

# load gpkgs
all_gpkg.ls <- lapply(all_gpkg, function(name) {vect(name)})


# merge them all in one file
merged_gpkg <- do.call("rbind", all_gpkg.ls)


crs(merged_gpkg)

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
  dplyr::mutate(ID = paste(country, region, group, point, sep="_")
  ) %>%  #dplyr::select(5:323) %>% names()
  dplyr::select(ID, region, country)

# convert back to terra object
xy_clean <- vect(data_clean)

nrow(xy_clean)
# loop over countries to run it faster -----------------------------------------

# select climate space variables: temp (t2m) and precipitation (tp)
clim_vars <- c("t2m", "tp")
temp_convert = 273.15  # convert temperature from Kelvin to Celsius

# Read .nc data as a raster in terra - way faster!  
dat_ras <- terra::rast("rawData/ERA_NET/data2.nc")


(var_names <- varnames(dat_ras))    # 7 vars: "d2m"   "t2m"   "tp"    "swvl1" "swvl2" "swvl3" "swvl4"
(ras_time  <- time(dat_ras))        #
(ras_source<- sources(dat_ras))        #
(n_layer   <- nlyr(dat_ras)) 

extract_clim_data <- function(country_name, ...) {
  print(paste("Processing", country_name))
  #
  #country_name = 'germany'
  
  #"11_25_116_4"
  
  #country <- vect(paste0('outData/dat_', country_name, '.gpkg'))
  
  #country_sf <- st_as_sf(country) 
  #country_sf %>% 
  #  filter(ID == "11_25_116_4")
  
  
  
  country <- vect(paste0('outData/dat_', country_name, '.gpkg'))
  xy_clean <- project(country, "EPSG:3035")
  
  # Get spatial data for each trap
  xy        <- terra::project(xy_clean, crs(dat_ras))
  
  # check projection
  crs(xy) == crs(dat_ras) 
  
  # get unique IDs and names
  xy$name = paste0(1:nrow(xy)) #1:nrow(xy_proj)
  
  # extract all values to the xy coordinates:
  dat_ext_df <- terra::extract(dat_ras, xy)
  
  # add group indication
  dat_ext_df$name <- xy$name
  dat_ext_df$ID   <- xy$ID
  
  # Create a datatable for each site (ID), variable, and time
  dat_ext_df <- as.data.table(dat_ext_df)
  
  # melt from wide to long format
  df_melt <- dat_ext_df %>%
    melt(id.vars = c('ID','name')) #', 

  # Add time  to df and split in months:
  df <- 
    df_melt %>% 
    arrange(ID, name, variable) %>%
    dplyr::mutate(time = rep(ras_time, nrow(dat_ext_df))) %>% 
    tidyr::separate(variable, 
                    c("var", "time_num", 'xx'), "_") %>% 
    na.omit() %>% 
    dplyr::filter(var %in% clim_vars ) %>% # filter only clim vars
    dplyr::mutate(year  = lubridate::year(time), 
                  month = lubridate::month(time), 
                  doy   =  lubridate::yday(time) + 1) # %>%  # as POXIT data has January 1st at 0
  
  # keep only means per year
  df_mean <- 
    df %>% 
    #na.omit(.) %>% 
    dplyr::select(-c(time_num, xx)) %>% # remove unnecessary cols
    spread(var, value) %>%
    mutate(t2m = t2m - temp_convert) %>% 
    group_by(ID, name, year) %>% 
    dplyr::summarize(tmp = mean(t2m),# mean annual temp
              prec = mean(tp)*1000*365)  # sum annual precip (it is mothly means), convert from meters to mm
  
  return(df_mean)
}



# run for all---------------------------------------------------------------------------------------------------

# list all countries
country_names <- list( "austria","belgium", "czechia", "france", "germany", "italy", "luxembourg", "poland", 
                       "slovakia", "slovenia", "switzerland") #austria" 



# loop over all countries
all_climate <- lapply(country_names, function(cn) {
  extract_clim_data(cn)
})


# Combine results from all countries
final_climate <- do.call(rbind, all_climate)

reference_period = 1980:2015


# Calculate temp anomalies
final_climate <- final_climate %>% 
  group_by(ID) %>%
  mutate(tmp_mean = mean(tmp[year %in% reference_period]),
         tmp_sd = sd(tmp[year %in% reference_period]),
         tmp_z  = (tmp - mean(tmp[year %in% reference_period])) / sd(tmp[year %in% reference_period]),
         prcp_z = (prec - mean(prec[year %in% reference_period])) / sd(prec[year %in% reference_period]))


length(unique(final_climate$ID))

final_climate_out <- final_climate %>% 
  ungroup(.) %>% 
  dplyr::select(-name) %>% 
  dplyr::filter(year %in% 2018:2020) %>% 
  dplyr::select(-tmp_mean, -tmp_sd)

length(unique(final_climate$ID))



# Export data -------------------------------------------------------------
fwrite(final_climate_out, 'outData/climate.csv')


save(
  final_climate,
  file="outData/plots_env.Rdata")













# test on single country -------------------------------------------------------
# Read field data
country_name = 'slovakia'

country <- vect(paste0('outData/dat_', country_name, '.gpkg'))
xy_clean <- project(country, "EPSG:3035")

# Read .nc data as a raster in terra - way faster!  
dat_ras <- terra::rast("rawData/ERA_NET/data.nc")

# Get spatial data for each trap
xy        <- terra::project(xy_clean, crs(dat_ras))

# check projection
crs(xy) == crs(dat_ras) 

# get unique IDs and names
xy$name = paste0(1:nrow(xy)) #1:nrow(xy_proj)

# extract all values to the xy coordinates:
dat_ext_df <- terra::extract(dat_ras, xy)

# add group indication
dat_ext_df$name <- xy$name
dat_ext_df$ID   <- xy$ID

# Create a datatable for each site (ID), variable, and time
dat_ext_df <- as.data.table(dat_ext_df)

# melt from wide to long format
df_melt <- dat_ext_df %>%
  melt(id.vars = c('ID','name')) #', 


# select climate space variables: temp (t2m) and precipitation (tp)
clim_vars <- c("t2m", "tp")
temp_convert = 273.15  # convert temperature from Kelvin to Celsius

# Add time  to df and split in months:
df <- 
  df_melt %>% 
  arrange(ID, name, variable) %>%
  dplyr::mutate(time = rep(ras_time, nrow(dat_ext_df))) %>% 
  tidyr::separate(variable, 
                  c("var", "time_num", 'xx'), "_") %>% 
  na.omit() %>% 
  dplyr::filter(var %in% clim_vars ) %>% # filter only clim vars
  dplyr::mutate(year  = lubridate::year(time), 
                month = lubridate::month(time), 
                doy   =  lubridate::yday(time) + 1) # %>%  # as POXIT data has January 1st at 0

# keep only means per year
df_mean <- 
  df %>% 
  na.omit(.) %>% 
  #dplyr::filter(var %in% clim_vars ) %>% # filter only clim vars
  dplyr::select(-c(time_num, xx)) %>% # remove unnecessary cols
  spread(var, value) %>%
  mutate(t2m = t2m - temp_convert) %>% 
  group_by(ID, name, year) %>% 
  summarize(tmp = mean(t2m),# mean annual temp
            prec = mean(tp)*1000*365)  # sum annual precip (it is mothly means), convert from meters to mm
















