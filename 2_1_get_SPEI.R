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
pattern_years = "^19(8[0-9]|9[0-9])|^20(0[0-9]|1[0-9]|2[0-1]).*\\.gz$"


# list all countries
country_names <- list( "austria","belgium", "czechia", "france", "germany", "italy", "luxembourg", "poland", 
                       "slovakia", "slovenia", "switzerland") #austria" 
# read global data
climate   <- fread("outData/climate_1980_2023.csv")


# test
# read local data
country_name = 'slovakia'
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

# get lat long values
crds(xy_latlng)

# get xy as wgs for teh SPEI
xy$x_wgs <- crds(xy_latlng)[,'x']
xy$y_wgs <- crds(xy_latlng)[,'y']


# add x_wgs to cluster data
df_clim1 <- df_clim %>% 
  left_join(as.data.frame(xy), by = join_by(ID)) %>% 
  dplyr::select(ID, year, tmp, prec, y_wgs)  

# !!! need to first get data by months and year!!!
  
head(df_clim1)


# get SPEI -----------------------------------------------------------
# SPEI: Standardized potential evapotranspiration index --------------------------------

# Function to calculate PET using Thornthwaite method
calculate_pet_thornthwaite <- function(temp, latitude) {
  # Using the built-in Thornthwaite function from the SPEI package
  PET <- thornthwaite(temp, latitude)
  return(PET)
}

#Get values for SPEI
df_th <-
  df_clim1 %>% 
  #na.omit() %>% # remove duplicated values
  group_by(ID) %>%
  mutate(PET = calculate_pet_thornthwaite(tmp, unique(y_wgs)[1]),
         BAL = prec - PET) 

# calculate SPEI for each falsto location:
df_ls <- df_th %>%
  #group_by(falsto_name, .add = TRUE) %>% 
  group_split(ID)


# Calculate the SPEI for each location:
get_SPEI <- function(df, ...){
  #df <- test
  # get XY name
  id = unique(df$ID)
  
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

