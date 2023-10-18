
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


#cz <- vect("C:/Users/ge45lep/Documents/2023_PanEuropean/r_paneurop/inData/data_CZ.gpkg")
cz2 <- vect("inData/data_CZ.gpkg")
plot(cz2)



