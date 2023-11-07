# get climatic characteristics
# temp
# prec
# temp anomalies???
# from ERA-NET data


library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(tidyr)
library(dplyr)
library(data.table)


library(PCICt)
library(zoo)                # for as.Date() specification

library(ggplot2)
library(ggpubr)

# get data from Franconia ------------------------------------------------
franconia <- vect("rawData/franconia/new_GPS/sites_final_updPassau.shp")

# get data from Pan-Europe 
#list all field data collected; get paths
all_gpkg <- list.files(path = "rawData/collected", 
                       pattern = "*.gpkg$", 
                       full.names = TRUE)  # read the full path to read the files properly

# load gpkgs
all_gpkg.ls <- lapply(all_gpkg, function(name) {vect(name)})

# merge them all in one file
merged_gpkg <- do.call("rbind", all_gpkg.ls)

crs(merged_gpkg)


# Load ERA5 data --------------------------------------------

# Data was downloaded as netCDF from https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land-monthly-means?tab=overview

# Read .nc data as a raster in terra - way faster!  
dat_ras <- terra::rast("rawData/ERA_NET/data.nc")

# Get spatial data for each trap
xy_proj        <- terra::project(merged_gpkg, crs(dat_ras))
franconia_proj <- terra::project(franconia, crs(dat_ras))
      
# correct                    
crs(xy_proj) == crs(dat_ras)
crs(franconia_proj) == crs(dat_ras)

# get unique IDs and names
xy_proj$ID = paste0('xy', 1:nrow(xy_proj)) #1:nrow(xy_proj)
#xy_proj$name = 'xy'

franconia_proj$ID = paste0('fr', 1:nrow(franconia_proj))
#franconia_proj$name = 'franconia'


# simlify spatial data: keep only location, name and id;
# remove unnecesary columns
xy    <- xy_proj[c("ID")] #
franc <- franconia_proj[c("ID")] # , "name"

# merge both shps together
xy <- rbind(xy, franc)

# extract all values to the xy coordinates:
dat_ext_df <- terra::extract(dat_ras, xy, df = TRUE)

# add group indication
dat_ext_df$name <- xy$ID
# Check the variables of interest:
(var_names <- varnames(dat_ras))    # 7 vars: "d2m"   "t2m"   "tp"    "swvl1" "swvl2" "swvl3" "swvl4"
(ras_time  <- time(dat_ras))        #
(ras_source<- sources(dat_ras))        #
(n_layer   <- nlyr(dat_ras))        #



# Create a datatable for each site (ID), variable, and time
dat_ext_df <- as.data.table(dat_ext_df)

# melt from wide to long format
df_melt <- dat_ext_df %>%
  melt(id.vars = c('ID','name')) #', 

#df_melt <- melt(setDT(dat_ext_df), id.vars = c('ID', 'name'))


# Add time  to df and split in months:
df <- 
  df_melt %>% 
  arrange(ID, variable) %>%
  dplyr::mutate(time = rep(ras_time, nrow(dat_ext_df))) %>% 
  tidyr::separate(variable, 
           c("var", "time_num", 'xx'), "_") %>% 
  na.omit() %>% 
  dplyr::mutate(year  = lubridate::year(time), 
                month = lubridate::month(time), 
               # day   = lubridate::day(time),
                doy   =  lubridate::yday(time) + 1) # %>%  # as POXIT data has January 1st at 0
  #dplyr::select(-day) 

# check values for precipitation?
df %>% 
  filter(var == 'tp'& ID == 5) %>%   # in meters depth, and day!!, monthly averages -> multiply the average by 365 to get yearly value! 
  na.omit() %>%                      # see days per each month for more precice: https://stackoverflow.com/questions/6243088/find-out-the-number-of-days-of-a-month-in-r 
  group_by(year) %>% 
  mutate(prcp = value*1000) %>% 
  summarize(sum_prcp = mean(prcp)*365)


# select climate space variables: temp (t2m) and precipitation (tp)
clim_vars <- c("t2m", "tp")
temp_convert = 273.15  # convert temperature from Kelvin to Celsius



# keep only means per year
df_mean <- 
  df %>% 
  na.omit(.) %>% 
  dplyr::filter(var %in% clim_vars ) %>% # filter only clim vars
  dplyr::select(-c(time_num, xx)) %>% # remove unnecessary cols
  spread(var, value) %>%
  mutate(t2m = t2m - temp_convert) %>% 
  group_by(ID, year) %>% 
  summarize(tmp = mean(t2m),# mean annual temp
            prec = mean(tp)*1000*365)  # sum annual precip (it is mothly means), convert from meters to mm
  

head(df_mean)

# classify warm years
df_mean2 <- df_mean %>% 
  mutate(cat = case_when(year %in% c(2018:2020) ~ 'dist',
                         .default = "ref")) #%>% 
  #group_by(ID, cat) %>% 
  #summarize(tmp = mean(tmp),# mean annual temp
   #         prec = sum(prec))  # sum annual precip, convert from meters to mm

  
df_mean_sub <- df_mean2 %>%
  ungroup(.) %>% 
  slice_sample(n = 1000)
  #dplyr::sample_n(2) #(n = 1000)# %>% # subset number of rows to speed up lotting
  

# get data fr reference
df_ref <- df_mean2 %>%
  filter(cat == 'ref') %>% 
  slice_sample(n = 5000)
# get data for warm years
df_warm <- df_mean2 %>%
  filter(cat == 'dist') 
#dplyr::sample_n(2) #(n = 1000)# %>% # subset number of rows to speed up lotting




windows()
df_ref %>%
  ggplot(aes(x = tmp,
             y = prec,
             color = cat)) +
  stat_density_2d(
    aes(fill = after_stat(density)), # Use density directly
    geom = "raster", 
    contour = FALSE
  ) +
  
  scale_fill_gradient2(
    low = "white", 
   # mid = "white", 
    high = "forestgreen", 
   # midpoint = 4e-04
  ) +
 geom_point(data = df_warm, aes(x = tmp,
                         y = prec), color = 'black', size = 0.5) +
  theme_minimal()
