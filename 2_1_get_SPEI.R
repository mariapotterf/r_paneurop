# Get SPEI
# calculate SPEI for each country and each location

# Read libs  --------------------------------------------------------------

library(SPEI)
library(sf)
library(dplyr)
library(data.table)
library(tidyr)
library(rgdal)
library(raster)
library(lubridate)
library(fasterize)
library(ggpubr)
library(terra)
library(R.utils)
library(stringr)

# SPEI scales:
my_spei_scales = c(3) #c(3,12) # 3,6,12
# filter through years: >1980 ------------------------
# 
#pattern_years = "^19(8[0-9]|9[0-9])|^20(0[0-9]|1[0-9]|2[0-3]).*\\.gz$" # filter data 1980-2023


# list all countries
country_names <- list( "austria","belgium", "czechia", "france", "germany", "italy", "luxembourg", "poland", 
                       "slovakia", "slovenia", "switzerland") #austria" 
# read global data
climate   <- fread("outData/climate_1980_2023_months.csv")


# vars

# SPEI: Standardized precipitation evapotranspiration index --------------------------------
SPEI_vars <- c("t2m", "tp")

temp_convert = 273.15  # convert temperature from Kelvin to Celsius
#number_days <- 365/12  # multiply teh prec by teh number of days in a month to correct precipitation



# Function to calculate PET using Thornthwaite method
calculate_pet_thornthwaite <- function(temp, latitude) {
  # Using the built-in Thornthwaite function from the SPEI package
  PET <- thornthwaite(temp, latitude)
  return(PET)
}

# Calculate the SPEI for each location:
get_SPEI <- function(df, ...){
  #df <- test
  # get XY name
  id = unique(df$ID)
  
  # convert df to time series
  df.ts <- df %>% 
    arrange(year, month) %>% 
    ts(df, start = c(1980, 01), end=c(2023,12), frequency=12) 
  
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
    mutate(ID = rep(id, nrow(out_scales)))
  # (out_scales)
  return(out_scales)
  
}



process_spei <- function(country_name, ...) {
  
  print(paste("Processing", country_name))
  # read local data
  #country_name = 'slovakia'
  # loop over names inside the function
  xy        <- vect(paste0('outData/dat_', country_name, '.gpkg'))
  
  
  # filter cluster for each country
  unique_cluster <- unique(xy$ID)
  
  
  df_clim <- climate %>% 
    filter(ID %in%unique_cluster)
  
  # Get spatial data for each plot
  # project to long lat for Thornweite for SPEI
  xy_latlng <- terra::project(xy, 
                              "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  # get xy as wgs for teh SPEI
  xy$y_wgs <- crds(xy_latlng)[,'y']
  
  
  #Get values for SPEI
  df <-
    df_clim %>% 
    dplyr::filter(var %in% SPEI_vars ) %>% # filter only soil water content
    na.omit() %>% # remove duplicated values
    dplyr::select(-c(time_num, xx)) %>% # remove unnecessary cols
    spread(var, value) %>%
    mutate(t2m = t2m - temp_convert,
           tp = tp*number_days)  %>% 
    rename(tmp = t2m,
           prcp = tp)
  
  # add x_wgs to cluster data
  df <- df %>% 
    left_join(as.data.frame(xy), by = join_by(ID)) %>% 
    dplyr::select(ID, year, month, tmp, prcp, y_wgs) %>% 
    arrange(year, month)
  
  
  
  
  # get thornthwaite -----------------------------------------------------------
  df_th <-
    df %>% 
    #na.omit() %>% # remove duplicated values
    group_by(ID) %>%
    mutate(PET = calculate_pet_thornthwaite(tmp, unique(y_wgs)[1]),
           BAL = prcp - PET) 
  
  # calculate SPEI for each location ------------------------------------------
  df_ls <- df_th %>%
    #group_by(falsto_name, .add = TRUE) %>% 
    group_split(ID)
  
  
  # apply over the list of df (locations):
  df_ls2<- lapply(df_ls, get_SPEI)
  
  # merge into one file:
  df_spei_ID <- do.call('rbind', df_ls2)
  
  return(df_spei_ID)
}


# loop over all countries
all_SPEI <- lapply(country_names, function(cn) {
  process_spei(cn)
})


# Combine results from all countries
final_SPEI <- do.call(rbind, all_SPEI)

length(unique(final_SPEI$ID))

fwrite(final_SPEI, 'outData/SPEI.csv')

