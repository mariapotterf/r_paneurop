# Extract future species distribution 
# wessely 2024

# process;

# read the list of all Wessely species (67)
# read my species list (~37)
# narrows down speies list
# read the vect points

# look over the rasters to extract the values
# this will be 0/1 per species and decade
# merge 

library("terra")
library(tidyr)
library(ggplot2)
library(viridis)
library(raster)
library(reshape2)



# Get species Wessely
species_wessely <- read.csv("rawData/wessely_species_bottleneck/data_submit/data_submit/species_list.csv")[,1]

