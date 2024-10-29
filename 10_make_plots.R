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

