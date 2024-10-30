# make plots:

# read final tables, models
# make fina plots and tables to export

gc()
# Libs --------------------------------------------------------------------------


library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

library(ggpubr)

# library
library(ggridges)

library(purrr)
library(scales)  # add % for the plots
library(RColorBrewer)


source('my_functions.R')

# Set a global theme -------------------------------------------------------------------

theme_set(
  theme_classic() + 
    theme(
      legend.position = 'bottom',
      text = element_text(size = 4),         # Set all text to size 8
      axis.text = element_text(size = 8),    # Axis tick labels
      axis.title = element_text(size = 8),   # Axis titles
      strip.text = element_text(size = 8),   # Facet labels
      legend.text = element_text(size = 8),  # Legend text
      plot.title = element_text(size = 8)    # Plot title
    )
)


# read data --------------------------------------------------------------------------

df_fin <- fread('outTable/df_fin.csv')  # update table with moe columsn 

# get vegetation data
load("outData/veg.Rdata")



# Get gables

#View(df_fin)
# get summary for the scatter plot, one point is one country & management
df_summary <- df_fin %>%
  # group_by(manag) %>% # country, 
  summarize(
    rich_med = median(rIVI, na.rm = TRUE),
    rich_sd = sd(rIVI, na.rm = TRUE),
    rich_25 = quantile(rIVI, 0.25, na.rm = TRUE), # rich_mean -rich_sd, #
    rich_75 = quantile(rIVI, 0.75, na.rm = TRUE), # rich_mean +rich_sd,
    dens_med = median(stem_density   , na.rm = TRUE),
    dens_sd   = sd(stem_density  , na.rm = TRUE),
    dens_25 = quantile(stem_density  , 0.25, na.rm = TRUE), # dens_mean - dens_sd, #
    dens_75 = quantile(stem_density  , 0.75, na.rm = TRUE), # dens_mean + dens_sd, #
    .groups = 'drop'
  )

(df_summary)

