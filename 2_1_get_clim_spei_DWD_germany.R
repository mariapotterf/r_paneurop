# Get climate DWD - only for germany - has better resolution 

# Get SPEI data

# at daily resolution
# for each trap

# SPEI data monitor the drought worldwide:

# Calculation of the Standardised Precipitation-Evapotranspiration Index
# Global SPEI database


# The Global SPEI database, SPEIbase, offers long-time, 
# robust information about drought conditions at the global scale, 
# with a 0.5 degrees spatial resolution and a monthly time resolution. 
# It has a multi-scale character, providing SPEI time-scales between 1 and 48 months. 
# Currently it covers the period between January 1901 and December 2020.

# Beguería S. (2017) SPEIbase: R code used in generating the SPEI global database, doi:10.5281/zenodo.834462.

gc()
rm(list = ls())

#spring.months         = 3:5
#veg.months            = 4:9 #3:10
#study.period.extended = 2012:2021




# Read my paths -----------------------------------------------------------
#source('myPaths.R')


# Read libs  --------------------------------------------------------------

library(SPEI)
library(sf)
library(dplyr)
library(data.table)
library(tidyr)
library(rgdal)
library(raster)
#library(tidyverse)
library(lubridate)
#library(patchwork)
library(fasterize)
library(ggpubr)
library(terra)
library(R.utils)
library(stringr)
library(terra)

# SPEI scales:
my_spei_scales = c(1, 6, 12) #c(3,12) # 3,6,12
# filter through years: >1980 ------------------------
# 
pattern_years = "^19(8[0-9]|9[0-9])|^20(0[0-9]|1[0-9]|2[0-3]).*\\.gz$"

# Get spatial data for each plot ---------
xy_subplots        <- vect("outData/dat_germany.gpkg") # read DE trap location


# converet to average coordinates per cluster ----------
# Example: If `xy` is a SpatVector with geometry
xy_sub_df <- as.data.frame(xy_subplots)  # Convert SpatVector to DataFrame

# Extract coordinates
coords <- crds(xy_sub_df)  # Extract X, Y coordinates
xy_sub_df$X <- coords[,1]
xy_sub_df$Y <- coords[,2]

# Compute mean coordinates per cluster
mean_coords <- xy_sub_df %>%
  group_by(cluster) %>%
  summarise(mean_X = mean(X, na.rm = TRUE),
            mean_Y = mean(Y, na.rm = TRUE))

# Convert back to SpatVector if needed
xy <- vect(mean_coords, geom = c("mean_X", "mean_Y"), crs = crs(xy))

# Print result
print(mean_coords)


# convert coordinate system
xy2       <- terra::project(xy, "EPSG:31467")  # coordinate system from the DWD data: Germany
xy_latlng <- terra::project(xy2, 
                            "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# get lat long values
crds(xy_latlng)

# keep only usnique localtions - 158 traps
xy2 <- xy2[, !names(xy2) %in% c('year', 'id')] # remove columns
xy2 <- unique(xy2)

# get xy as wgs for teh SPEI
xy2$x_wgs <- crds(xy_latlng)[,'x']
xy2$y_wgs <- crds(xy_latlng)[,'y']


# stack rasters firsst
process_raster_file <- function(file, path_prefix, var) {
  ras_path <- paste(path_prefix, var, file, sep = "/")
  unzipped_path <- gsub(".gz$", "", ras_path)
  
  # Check if unzipped file already exists
  if (!file.exists(unzipped_path)) {
    R.utils::gunzip(ras_path, remove = FALSE, overwrite = TRUE)
  }
  
  z <- rast(unzipped_path)
  crs(z) <- "EPSG:31467"
  return(z)
}


result_list <- list()
vars <- c("temp", "precip")  # Assuming these are your variables

for (i in vars){
  file_ls <- list.files(paste("rawData", i, sep = "/"),
                        pattern = pattern_years,#"^19(8[0-9]|9[0-9])|^20(0[0-9]|1[0-9]|2[01]).*\\.gz$",
                        recursive=TRUE)
  
  #(file_ls)
  # additional filtering to keep only files with '.gz' at the end
  file_ls_gz_only <- file_ls[grepl("\\.gz$", file_ls)]
  #(file_ls_gz_only)
  
  # Assuming file_ls contains the paths of all raster files
  ras_ls <- lapply(file_ls_gz_only, process_raster_file, path_prefix = paste("rawData", sep = "/"), var = i)
  
  ras_stack <- rast(ras_ls)
  
  # extract values
  extracted_values <- terra::extract(ras_stack, xy2)
  
  # add plot ID
  extracted_values$cluster <- xy2$cluster
 # extracted_values$globalid    <- xy2$globalid
  
 # names(extracted_values) <- gsub('asc', '', names(extracted_values))
  
  long.df <- extracted_values %>%
    pivot_longer(!c(ID,cluster), names_to = "time", values_to = 'vals') %>%
    mutate(month = as.integer(str_sub(time, -2, -1)),
           year = as.integer(str_sub(time, 5,6))) %>%
    dplyr::select(c(time))
  
  result_list[[i]] <- long.df
  
}

# merge data together to calculate SPEI -----------------------------
df_prec <- result_list$precip
df_temp <- result_list$temp

df_prec <- df_prec %>% 
  rename(prcp = vals)

df_temp <- df_temp %>% 
  rename(tmp = vals) %>% 
  mutate(tmp = tmp/10)  # correct for DWD recording: values are in 1/10 C format


hist(df_temp$tmp)
hist(df_prec$prcp)


# merge data
df_clim <- df_prec %>% 
  left_join(df_temp, by = join_by(globalid, falsto_name, month, year, ID)) %>% 
  right_join(as.data.frame(xy2  ), by = c('falsto_name',  'globalid'                )) %>% 
  dplyr::select(c(falsto_name, prcp, month,  year,   tmp, y_wgs)) %>% 
  distinct()

summary(df_clim)
(df_clim)

unique(df_clim$falsto_name)


# get SPEI -----------------------------------------------------------
# SPEI: Standardized potential evapotranspiration index --------------------------------

# Function to calculate PET using Thornthwaite method
calculate_pet_thornthwaite <- function(temp, latitude) {
  # Using the built-in Thornthwaite function from the SPEI package
  PET <- thornthwaite(temp, latitude)
  return(PET)
}

#Get values for SPEI
df1 <-
  df_clim %>% 
  na.omit() %>% # remove duplicated values
  group_by(falsto_name) %>%
  mutate(PET = calculate_pet_thornthwaite(tmp, unique(y_wgs)[1]),
         BAL = prcp - PET) 

# calculate SPEI for each falsto location:
df_ls <- df1 %>%
  #group_by(falsto_name, .add = TRUE) %>% 
  group_split(falsto_name)

#df.ts$SPEI <- dd




# Calculate the SPEI for each location:
get_SPEI <- function(df, ...){
  #df <- test
  # get XY name
  id = unique(df$falsto_name)
  
  # convert df to time series
  df.ts <- df %>% 
    arrange(year, month) %>% 
    ts(df, start = c(1980, 01), end=c(2021,12), frequency=12) 
  #ts(df, start = c(2008, 01), end=c(2021,12), frequency=12) 
  
  # Calculate spei or different time intervals:
  spei_ls <- lapply(my_spei_scales, function(s) {
    
    # extract just values from SPEI object:
    dd = spei(df.ts[,'BAL'], scale = s)$fitted
    
    # covert to dataframe, convert date to format
    df.out <- data.frame(spei=as.matrix(dd), 
                         date=zoo::as.Date(time(dd)))
    
    # add scale indication
    df.out <-df.out %>% 
      mutate(scale = rep(s, nrow(df.out)))
    return(df.out)
  })
  # merge scaes tables
  out_scales = do.call('rbind', spei_ls)
  
  # add location indication
  out_scales <-out_scales %>% 
    mutate(falsto_name = rep(id, nrow(out_scales)))
  # (out_scales)
  return(out_scales)
  
}

# apply over the list of df (locations):
df_ls2<- lapply(df_ls, get_SPEI)

# merge into one file:
df_spei_ID <- do.call('rbind', df_ls2)

df_spei_ID2 <- df_spei_ID

df_spei_ID2 %>% 
  # filter(falsto_name == 'Bleichach_1') %>% 
  filter(scale == 12) %>%
  ggplot(aes(x = date, y = spei)) +
  #geom_line() +
  stat_summary(fun = mean, geom = "point", color = "black", size = 1) +
  scale_x_date(limits = as.Date(c("1980-01-01", "2021-12-01"))) # Replace with actual dates


hist(df_prec$prcp)


df_spei_ID2 %>% 
  filter(falsto_name == 'Blaichach_1') %>% 
  filter(scale == 12) %>%
  ggplot(aes(x = date, y = spei)) +
  geom_line() 


df_prec %>% 
  group_by(year, falsto_name) %>% 
  summarize(sum_prcp = sum(prcp)) %>%
  View()
ggplot(aes(x = year,
           y = sum_prcp)) +
  geom_line(alpha = 0.5)

# summarize spei per year - one SPEI value per year and ID! 
df_spei_ID <-  df_spei_ID %>% 
  dplyr::mutate(year  = lubridate::year(date),
                month = lubridate::month(date)) %>%  
  dplyr::select(-c(date))


# plot if correct??
df_spei_ID %>%
  filter(year %in% 2015:2021) %>% 
  ggplot(aes(x = month ,
             y = spei,
             group = factor(year),
             color = factor(year))) +
  
  #geom_point(alpha = 0.5)  +
  geom_smooth()+
  #geom_point(alpha = 0.5)  +
  facet_grid(year~scale) +
  geom_hline(yintercept = 0, col = 'red', lty = 'dashed') 





# export file:
#fwrite(out.df, paste(myPath, outTable, 'xy_spei.csv', sep = "/"))

# spei does only up to 2021
df_spei_veg_season <- 
  df_spei_ID %>% 
  dplyr::filter(month %in% veg.months & year %in% study.period.extended) %>% 
  ungroup(.) %>% 
  group_by(falsto_name, year, scale) %>% 
  summarise(spei = median(spei)) 




# Export data -------------------------------------------------------------


data.table::fwrite(df_clim, 
                   'outTable/xy_clim_DWD.csv')
data.table::fwrite(df_spei_veg_season, 
                   'outTable/xy_spei_veg_season_DWD.csv')

data.table::fwrite(df_spei_ID, 
                   'outTable/xy_spei_all_DWD.csv')

