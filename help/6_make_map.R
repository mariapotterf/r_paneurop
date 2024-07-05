# Make map

# Load the packages
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(sf)
library(terra)

dat <- vect('rawData/outField_noItaly.gpkg')
hotspots <- vect('rawData/EU_hotspots.shp')


# Convert to sf objects for compatibility with ggplot2
dat_sf <- st_as_sf(dat)
hotspots_sf <- st_as_sf(hotspots)

# Get world map data and filter for Europe
world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- world[world$continent == "Europe",]

# List of European countries to highlight
highlight_countries <- c("France", "Germany", "Belgium", "Czechia", 
                         "Slovakia", "Slovenia", "Austria", "Poland", "Luxembourg")
europe$highlight <- ifelse(europe$name %in% highlight_countries, "Highlighted", "Not Highlighted")

# Define the coordinate reference system (CRS) for the map
crs_europe <- st_crs(3035) # EPSG:3035

# Transform all spatial data to the new CRS
europe <- st_transform(europe, crs_europe)
dat_sf <- st_transform(dat_sf, crs_europe)
hotspots_sf <- st_transform(hotspots_sf, crs_europe)


crs(europe)
crs(dat_sf)
crs(hotspots_sf)


# Plot the map with additional layers
ggplot() +
  geom_sf(data = europe, aes(fill = highlight), color = "white") +
  geom_sf(data = dat_sf, color = "red", size = 2) + # Adjust color and size as needed
  geom_sf(data = hotspots_sf, color = "green", size = 2) + # Adjust color and size as needed
  scale_fill_manual(values = c("Highlighted" = "blue", "Not Highlighted" = "lightgray")) +
  coord_sf(crs = crs_europe, xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal() +
  labs(title = "Map of Europe with Selected Countries and Additional Layers", fill = "Country Status")
