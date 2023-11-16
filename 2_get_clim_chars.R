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
library(stringr)

# get data from Franconia ------------------------------------------------
franconia <- vect("rawData/franconia/new_GPS/sites_final_updPassau.shp")

franconia_clim_space <- vect("rawData/franconia/extra/climate_hotspots2.gpkg")

# get data from Pan-Europe 
#list all field data collected; get paths
all_gpkg <- list.files(path = "rawData/collected", 
                       pattern = "*.gpkg$", 
                       full.names = TRUE)  # read the full path to read the files properly


# remove Italy, Szitzerland, slovenia
filtered_gpkg <- all_gpkg[!grepl("it_|sw_", all_gpkg)]

# load gpkgs
all_gpkg.ls <- lapply(filtered_gpkg, function(name) {vect(name)})

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
franconia_proj <- franconia_proj[order(franconia_proj$Name), ]

franconia_clim_space_proj <- terra::project(franconia_clim_space, crs(dat_ras))
# correct                    
crs(xy_proj) == crs(dat_ras)
crs(franconia_proj) == crs(dat_ras)
crs(franconia_clim_space_proj) == crs(dat_ras)

# get unique IDs and names
xy_proj$name = paste0('xy', 1:nrow(xy_proj)) #1:nrow(xy_proj)
franconia_proj$name = paste0('fr', 1:nrow(franconia_proj))
franconia_clim_space_proj$name = paste0('extra', 1:nrow(franconia_clim_space_proj))



# simlify spatial data: keep only location, name and id;
# remove unnecesary columns
xy    <- xy_proj[c("name")] #
franc <- franconia_proj[c("name")] # , "name"
franc_extra <- franconia_clim_space_proj[c("name")] # , "name"

# merge both shps together
xy <- rbind(xy, franc, franc_extra)

# extract all values to the xy coordinates:
dat_ext_df <- terra::extract(dat_ras, xy)

# add group indication
dat_ext_df$name <- xy$name
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
  arrange(ID, name, variable) %>%
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
  filter(var == 'tp'& name == 'xy5') %>%   # in meters depth, and day!!, monthly averages -> multiply the average by 365 to get yearly value! 
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
  group_by(ID, name, year) %>% 
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


# get data fr reference
df_ref <- df_mean2 %>%
  filter(cat == 'ref') %>% 
  filter(str_starts(name, "xy")) #%>% 


# get data for warm years
df_dist <- df_mean2 %>%
  filter(year %in% c(2018:2020))

# get data for warm years !!!!
df_franc <- df_mean2 %>%
  filter(cat == 'dist')  %>% 
  filter(str_starts(name, "fr")) %>% 
  ungroup() %>% 
  mutate(grp = factor(rep(1:42, each = 9)))



windows()

rescale_density <- function(x) {
  scales::rescale(x, to = c(0.001, 1))
}

# get summary table for franconia
df_summary <- df_franc %>%
  group_by(grp) %>%
  summarize(
    tmp_mean = mean(tmp),
    tmp_min = min(tmp),
    tmp_max = max(tmp),
    prec_mean = mean(prec),
    prec_min = min(prec),
    prec_max = max(prec),
    .groups = 'drop'
  )


# get final plot --------------------------------------------------------------

library(MASS)

# Function to calculate levels
calculate_levels <- function(density, percentages) {
  sorted_density <- sort(density, decreasing = TRUE)
  cum_density <- cumsum(sorted_density) / sum(sorted_density)
  sapply(percentages, function(p) {
    sorted_density[max(which(cum_density <= p))]
  })
}

# Calculate 2D density
d <- kde2d(df_dist$tmp, df_dist$prec, n = 500)
density_values <- d$z

# Set densities outside the range of your data to zero
density_values[density_values < 1e-6] <- 0

# Calculate the levels for specified percentages
levels <- calculate_levels(as.vector(density_values), c(0.50, 0.75, 0.90, 0.99, 1))

# Prepare data for ggplot
plot_data <- expand.grid(tmp = d$x, prec = d$y)
plot_data$density <- as.vector(density_values)

# Use cut to create factor levels, including one for zero density
plot_data$level <- cut(plot_data$density, breaks = c(-Inf, levels, Inf), labels = FALSE, include.lowest = TRUE)

# Define colors (from red to yellow, plus white for zero density)

library(RColorBrewer)
blue_colors <- brewer.pal(5, "Greys")

# Add 'white' at the beginning
my_colors <- c("white", blue_colors)


# Create the plot
p<-
  ggplot(plot_data) +
  geom_raster(aes(x = tmp, y = prec, fill = factor(level)), alpha = 0.6) +
  scale_fill_manual(values = my_colors,
                    labels = c("", rev(c("50%", "75%", "90%", "99%", "100%"))),
                    name = "Density") +
  geom_point(data = df_summary, 
             aes(x = tmp_mean, y = prec_mean, group = grp), color = 'white', size = 2) +
  geom_point(data = df_summary, 
             aes(x = tmp_mean, y = prec_mean, group = grp), color = 'black', size = 1.5) +
  geom_errorbar(data = df_summary, 
                aes(x = tmp_mean, 
                    #y = prec_mean, 
                    ymin = prec_min, 
                    ymax = prec_max, group = grp), 
                width = 0.1, 
                linewidth = 0.2,
                color = 'black') +
  geom_errorbarh(data = df_summary, 
                 aes(#x = tmp_mean, 
                     y = prec_mean, 
                     xmin = tmp_min, 
                     xmax = tmp_max, group = grp),
                 height =20, 
                 linewidth = 0.2,
                 color = 'black') +
  ylab('Mean annual precipitation [mm]') +
  xlab(bquote('Mean annual temperature [' * degree * 'C]')) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_rect(fill = "white", colour = NA),
        aspect.ratio = 1,
        axis.title.x = element_text(size = 10),  # X-axis label
        axis.title.y = element_text(size = 10),  # Y-axis label
        legend.title = element_blank(),  # Legend title
        legend.text = element_text(size = 10),   #   # Legend item text)
        legend.key.size = unit(0.3, "cm"))  

(p)

ggsave("franc_in_clim_space.pdf", 
       plot = p, 
       device = "pdf", 
       units = "cm",
       dpi = 300,
       width = 12, height = 10)

ggsave("franc_in_clim_space.png", 
       plot = p, 
       device = "png", 
       units = "cm",
       dpi = 300,
       width = 12, height = 10)

ggsave("franc_in_clim_space.svg", 
       plot = p, 
       device = "svg", 
       units = "cm",
       dpi = 300,
       width = 10, height = 8)


