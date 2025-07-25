# species climate suitability

# Extract future species distribution from Wessely 2024

# process;

# read the list of all Wessely species (67)
# read my species list (~37)
# narrows down species list
# read the vect points

# look over the rasters to extract the values
# this will be 0/1 per species and decade
# merge 

rm(list = ls())             # Remove all objects from global environment
graphics.off()              # Close all open graphics devices
gc()  

# Libs --------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
library(RColorBrewer)
library(cowplot)
library(colorspace)
library(here)
library(ggalluvial)
library(sf)
library(sjPlot)



# Input data -------------------------------------------------------------

source(file.path(public_dir, "code", "00_paths_functions.R"))

## Set a global theme -------------------------------------------------------------------

theme_set(
  theme_classic() + 
    theme(
      legend.position = 'bottom',
      text = element_text(size = 4),         # Set all text to size 8
      axis.text = element_text(size = 7),    # Axis tick labels
      axis.title = element_text(size = 7),   # Axis titles
      strip.text = element_text(size = 7),   # Facet labels
      legend.text = element_text(size =7),  # Legend text
      plot.title = element_text(size = 7)    # Plot title
    )
)





# 7. Species climate suitability  -------------------------------------------

df_stem_species_class %>% 
  dplyr::filter(stem_density > 0) %>% 
  distinct(Species) %>% 
  count()  #@ 35

# compare current species composition (from field data) with climate suitability from wessely
# Need to adjust tree species naming a bit: 
#   we have less species then wessely: 
# - eg we have besp, which can be betula pendula, betula pubescense
# - same for populus (3), quercus (6), salix ...


# identify what species are present per plot - use simple look up table, 
# add only latin names
species_present_field <- 
  df_stem_species_class %>% 
  dplyr::filter(VegType != "Mature") %>% 
  ungroup() %>% 
  group_by(plot, Species) %>% 
  summarize(reg_stem_density = sum(stem_density, na.rm = T )) %>% 
  mutate(presence = if_else(reg_stem_density > 0, 1, 0)) %>% 
  # dplyr::select(-reg_stem_density ) %>% 
  left_join(lookup_species_acronyms, by = c("Species" = "acc")) %>% 
  dplyr::select(-latin)  %>% 
  dplyr::rename(acc = Species) 


# add acronyms and consider presence across several species: eg betula sp.
species_suitability_summary <- 
  species_climate_suitability_raw  %>%
  rename(plot = site) %>% 
  # add naming
  left_join(lookup_species_field_wessely, by = c('species' = 'wessely')) %>%
  # group by species to allow occurence of species that have specified genus: eg betula sp.
  group_by(plot, scenario, acc) %>% 
  # Summarize and set sum_presence to 1 if the sum is greater than 1 - 
  # this account for the fact that wessely can have two betulas, I have only 1
  summarize(sum_presence = pmin(sum(overall_presence), 1), .groups = 'drop')

species_suitability_wide <- species_suitability_summary %>%
  pivot_wider(names_from = scenario, values_from = sum_presence) 

# merge both tables: the presently recorded species and species under climate scenarios
df_suitability_observed_vs_modeled <- species_suitability_wide %>% 
  left_join(species_present_field) %>%   # use left join to explude species that are recorded in field, but not present in Wessely database
  dplyr::rename(current = presence)


## Combine field occurence with vs Wessely species databases  ---------
length(unique(df_suitability_observed_vs_modeled$acc))  # final length is 30 species: corss between wessely and observed field database

# create look up table to identify the countries based on regions
region_country_lookup <- df_fin %>%
  group_by(country_pooled) %>%
  dplyr::summarise(region = list(unique(region))) %>%
  unnest(region) %>%
  mutate(region = as.integer(as.character(region)))

# add country indication to merged field and modeled suitability
df_suitability_observed_vs_modeled_reg <- 
  df_suitability_observed_vs_modeled %>%
  # Extract the first two characters of 'plot' as 'region' and convert to integer
  mutate(region = as.integer(substr(plot, 1, 2))) %>%
  # Left join with unique_regions_per_country to get country indication
  left_join(region_country_lookup)



# classify plot level species suitability 
df_suitability_observed_vs_modeled_reg <- df_suitability_observed_vs_modeled_reg %>%
  mutate(
    suitability_rcp26 = case_when(
      current == 1 & rcp26 == 1 ~ "suitable",
      current == 1 & rcp26 == 0 ~ "not_suitable",
      current == 0 & rcp26 == 1 ~ "novel",
      current == 0 & rcp26 == 0 ~ "0"
    ),
    suitability_rcp45 = case_when(
      current == 1 & rcp45 == 1 ~ "suitable",
      current == 1 & rcp45 == 0 ~ "not_suitable",
      current == 0 & rcp45 == 1 ~ "novel",
      current == 0 & rcp45 == 0 ~ "0"
    ),
    suitability_rcp85 = case_when(
      current == 1 & rcp85 == 1 ~ "suitable",
      current == 1 & rcp85 == 0 ~ "not_suitable",
      current == 0 & rcp85 == 1 ~ "novel",
      current == 0 & rcp85 == 0 ~ "0"
    )
  )

###  Overall analysis --------------------------
#### Share of stems that is not suitable under CC scenarios? --------------------------
# Convert suitability columns to logical: TRUE if unsuitable ("0"), FALSE otherwise

anyNA(df_suitability_observed_vs_modeled)  

# get total sum of stems: acdjust for fact that ionly 30 species overlaps between 2 dababases
total_stems         <- sum(df_stem_species_class$stem_density)
total_stems_wessely <- sum(df_suitability_observed_vs_modeled$reg_stem_density )
total_stems_wessely
total_stems
#6172500 - all stems, 5908250 over crossed databases 

# Convert to long format 
df_climate_suitability_regen_species_long <- df_suitability_observed_vs_modeled_reg %>%
  dplyr::select(plot, 
                acc,
                reg_stem_density , 
                country_pooled, 
                suitability_rcp26, 
                suitability_rcp45, 
                suitability_rcp85) %>%
  pivot_longer(cols = starts_with("suitability_rcp"), 
               names_to = "scenario", 
               values_to = "suitability") %>%
  mutate(scenario = gsub("suitability_", "", scenario)) #


df_climate_suitability_regen_species_long %>%  # Remove "suitability_" prefix
  group_by(suitability, scenario) %>% 
  dplyr::filter(suitability == "not_suitable") %>%  
  summarize(sum = sum(reg_stem_density)) %>% 
  mutate(share = round(sum/total_stems*100,1))

# suitability  scenario     sum share
# <chr>        <chr>      <int> <dbl>
#   1 not_suitable rcp26    3782000  61.3
# 2 not_suitable rcp45    4438000  71.9
# 3 not_suitable rcp85    5092750  82.5


#### Most affected species  ---------------
# # Identify species occurence per splots, and which species wuill lose teh most of climate suitability

species_plot_suitability_summary <- 
  df_climate_suitability_regen_species_long %>%
  dplyr::filter(reg_stem_density > 0) %>%
  
  # Count total number of unique plots per species
  group_by(acc) %>%
  mutate(n_plots_current_occurrence  = n_distinct(plot)) %>%
  # View()
  
  # Filter to only suitable cases
  dplyr::filter(suitability == "suitable") %>%
  
  # Count number of suitable plots per species and scenario
  distinct(plot, scenario, acc, n_plots_current_occurrence ) %>%
  count(acc, scenario, name = "n_suitable_plots") %>%
  
  # Join back total plot counts
  left_join(
    df_climate_suitability_regen_species_long %>%
      dplyr::filter(reg_stem_density > 0) %>%
      group_by(acc) %>%
      summarise(n_plots_current_occurrence = n_distinct(plot), .groups = "drop"),
    by = "acc"
  ) %>%
  
  # Calculate share
  mutate(n_unsuitable_plots = n_plots_current_occurrence - n_suitable_plots,
         share_suitable     = round(n_suitable_plots / n_plots_current_occurrence * 100,1),
         share_unsuitable   = 100 - share_suitable, #round(n_unsuitable_plots / n_plots_current_occurrence * 100,1),
         share_overall      = round(n_plots_current_occurrence / n_total_plots * 100,1)) %>%
  #View()
  dplyr::select(scenario, 
                share_overall,
                n_suitable_plots, 
                n_unsuitable_plots,     
                n_plots_current_occurrence,   # number of plots that species occured on
                share_suitable,         
                share_unsuitable,
                acc)

species_plot_suitability_summary

# identify the most affected species (highest share of unsuitable on plot in RCP45)
species_plot_suitability_summary %>%
  dplyr::filter(scenario == "rcp45") %>%
  #dplyr::filter(acc %in% top_species_plot_share) %>% 
  mutate(share_unsuitable = n_unsuitable_plots / n_plots_with_species * 100) %>%
  arrange(desc(share_overall)) %>% # , share_unsuitable
  View()




### Sankey plot test: need to restructure data -----------------

# update naming and coloring
species_colors_grey <- c(
  species_colors,
  other = "#D3D3D3"  # light grey
)

species_labels_grey <- c(
  species_labels,
  other = "Other species"
)


df_top_species_ovr <- 
  df_suitability_observed_vs_modeled_reg %>% 
  dplyr::mutate(
    acc = ifelse(acc %in% v_top10_species, acc, "other")  # relabel all others
  ) %>%
  #dplyr::filter(acc %in% top_species_plot_share) %>% 
  dplyr::filter(current == 1) %>% 
  dplyr::select(
    plot,
    acc,
    current,
    rcp26,
    rcp45, 
    rcp85,
    reg_stem_density,
    country_pooled
  ) %>% 
  group_by(acc) %>% 
  summarize( n_plots_present = n_distinct(plot),
             avg_stem_dens_current = sum(reg_stem_density)/n_total_plots,
             avg_stem_dens_rcp26 = sum(reg_stem_density[rcp26 == 1]) / n_total_plots,
             avg_stem_dens_rcp45 = sum(reg_stem_density[rcp45 == 1]) / n_total_plots,
             avg_stem_dens_rcp85 = sum(reg_stem_density[rcp85 == 1]) / n_total_plots,
             .groups = "drop") #%>% 


# Prepare long-format data
df_alluvial_stem_density_overall <- 
  df_top_species_ovr %>%
  dplyr::select(acc, 
                avg_stem_dens_current, 
                avg_stem_dens_rcp26, 
                avg_stem_dens_rcp45,
                avg_stem_dens_rcp85) %>%
  pivot_longer(
    cols = starts_with("avg_stem_dens"),
    names_to = "scenario",
    values_to = "stem_density"
  ) %>%
  mutate(
    scenario = case_when(
      scenario == "avg_stem_dens_current" ~ "Current",
      scenario == "avg_stem_dens_rcp26"   ~ "RCP2.6",
      scenario == "avg_stem_dens_rcp45"   ~ "RCP4.5",
      scenario == "avg_stem_dens_rcp85"   ~ "RCP8.5",
      TRUE ~ scenario
    )
  ) %>% 
  mutate(acc = factor(acc, levels = c(v_top10_species, 'other'))  # Set desired order
  )



# Create the alluvial plot
p_species_stem_density_ovr <- 
  ggplot(df_alluvial_stem_density_overall,
         aes(x = scenario,
             y = stem_density/1000,
             stratum = acc,
             alluvium = acc,
             fill = acc)) +
  geom_flow(alpha = 0.85) +
  geom_stratum() +
  theme_classic2(base_size = 8) +
  labs(title = "",
       x = "Scenario",
       #y =   y_lab_reg, #expression("Reg. stem density [n  1000 " " "*ha^{-1}*"]"),
       y = expression("Regeneration stem density [1000 n "*ha^{-1}*"]"),
       fill = "Species") +
  scale_fill_manual(values = species_colors_grey,
                    labels = species_labels_grey) +
  theme(
    legend.position = "right",
    legend.text = element_text(face = "italic", size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )

p_species_stem_density_ovr


ggsave(
  filename = file.path(public_dir, "figs", "Fig4.png"),
  plot = p_species_stem_density_ovr,
  width = 5,
  height = 3,
  dpi = 300,
  units = "in",
  bg = "white"  # Add white background to ensure clean rendering
)




### Get Tables:  -----------------------------------------------------------
# 1. share of stems/country
# 2. share of plots whene no stemps will remain
# merge into single table

# Summarise stem counts per plot, scenario, country
df_stem_suitability <- 
  df_climate_suitability_regen_species_long %>%
  group_by(country_pooled, scenario) %>%
  summarise(n_plot_country = n_distinct(plot),
            total_stems = sum(reg_stem_density, na.rm = TRUE),
            suitable_stems = sum(reg_stem_density[suitability == "suitable"], na.rm = TRUE),
            .groups = "drop"
  ) %>%
  mutate(
    share_suitable = ifelse(total_stems > 0, suitable_stems / total_stems, 0)  # <- Treat 0-stem plots as 0%
  )

# Step 2: Average share across plots (including zero-stem ones as 0%)
df_stem_suitability_avg <- df_stem_suitability %>%
  group_by(country_pooled, scenario, n_plot_country) %>%
  summarise(
    avg_share_suitable = round(mean(share_suitable) * 100, 1),  # As percentage
    .groups = "drop"
  ) %>% 
  mutate(
    country_full = case_when(
      country_pooled == "AT" ~ "Austria",
      country_pooled == "CH" ~ "Switzerland",
      country_pooled == "CZ" ~ "Czech Republic",
      country_pooled == "DE" ~ "Germany",
      country_pooled == "FR" ~ "France",
      country_pooled == "PL" ~ "Poland",
      country_pooled == "SI" ~ "Slovenia",
      country_pooled == "SK" ~ "Slovakia",
      TRUE ~ NA_character_  # Catch anything unexpected
    )
  ) %>%
  arrange(country_full)


df_stem_suitability_avg_wide <- df_stem_suitability_avg %>%
  select(country_pooled, country_full, scenario, avg_share_suitable) %>%
  pivot_wider(
    names_from = scenario,
    values_from = avg_share_suitable,
    names_prefix = "stem_share_"
  ) %>%
  arrange(country_full)


# Result
df_stem_suitability_avg_wide



#### Summary table per country -----------------------------

# get tables: 
# species richness per country 
#### Prepare Table 1: share of climatically suitable tree species per plot/country ------------------

# create a df that contains all plots (also empty ones to calculate properly averages)
df_master <- df_fin %>% 
  dplyr::select(plot, 
                country_pooled )


# compare teh shares of species that are climatically suitable per plot
df_suitability_richness <-   df_suitability_observed_vs_modeled_reg %>% 
  
  group_by(plot, country_pooled) %>%
  mutate(
    
    # calculate number of species
    richness_current = length(unique(acc[current == 1])),
    richness26       = length(unique(acc[current == 1 & rcp26 == 1])),
    richness45       = length(unique(acc[current == 1 & rcp45 == 1])),
    richness85       = length(unique(acc[current == 1 & rcp85 == 1])),
    
    # get shares                     
    rcp26_share = length(unique(acc[current == 1 & rcp26 == 1])) /
      length(unique(acc[current == 1])) * 100,
    rcp45_share = length(unique(acc[current == 1 & rcp45 == 1])) /
      length(unique(acc[current == 1])) * 100,
    rcp85_share = length(unique(acc[current == 1 & rcp85 == 1])) /
      length(unique(acc[current == 1])) * 100
  ) %>%
  ungroup() %>% 
  dplyr::select(plot,   
                acc,
                richness_current,
                richness26,
                richness45,
                richness85,
                current, 
                rcp26, rcp26_share, 
                rcp45, rcp45_share, 
                rcp85, rcp85_share ) %>%
  # filter only present species
  dplyr::filter(current == 1) %>%
  #  View()
  # merge back full table with empty plots
  full_join(df_master) %>% 
  replace_na(list(
    current = 0,
    rcp26_share = 0,
    rcp45_share = 0,
    rcp85_share = 0,
    richness_current = 0,
    richness26 = 0,
    richness45 = 0,
    richness85 = 0
  ))# %>%


#  Calculate shares 
df_suitability_richness_share <- df_suitability_richness %>% #  View()
  mutate(country_pooled = ifelse(is.na(country_pooled), "FR", 
                                 as.character(country_pooled))) %>%
  mutate(country_pooled = as.factor(country_pooled))  %>% # Optional: re-convert to factor if needed
  group_by(country_pooled) %>% 
  #  head()
  summarize(richness   = mean(richness_current),
            richness26 = mean(richness26),
            richness45 = mean(richness45),
            richness85 = mean(richness85),
            richness26_share = mean(rcp26_share),
            richness45_share = mean(rcp45_share),
            richness85_share = mean(rcp85_share)) %>% 
  mutate(
    country_full = case_when(
      country_pooled == "AT" ~ "Austria",
      country_pooled == "CH" ~ "Switzerland",
      country_pooled == "CZ" ~ "Czech Republic",
      country_pooled == "DE" ~ "Germany",
      country_pooled == "FR" ~ "France",
      country_pooled == "PL" ~ "Poland",
      country_pooled == "SI" ~ "Slovenia",
      country_pooled == "SK" ~ "Slovakia",
      TRUE ~ NA_character_  # Catch anything unexpected
    )
  ) %>%
  arrange(country_full)



#### Lost plots:  how many plots per country does not contain any currently present species??? ---------------

df_suitability_plots <- df_suitability_observed_vs_modeled_reg %>%
  group_by(country_pooled, plot) %>%
  mutate(
    any_present = any(current == 1),
    lost26 = if (any_present[1]) all(rcp26[current == 1] == 0) else FALSE,
    lost45 = if (any_present[1]) all(rcp45[current == 1] == 0) else FALSE,
    lost85 = if (any_present[1]) all(rcp85[current == 1] == 0) else FALSE
  ) %>%
  ungroup() %>%
  group_by(country_pooled) %>%
  mutate(
    n_plots_total = n_distinct(plot),
    n_lost_26 = n_distinct(plot[lost26]),
    n_lost_45 = n_distinct(plot[lost45]),
    n_lost_85 = n_distinct(plot[lost85]),
    share_lost_26 = n_lost_26 / n_plots_total * 100,
    share_lost_45 = n_lost_45 / n_plots_total * 100,
    share_lost_85 = n_lost_85 / n_plots_total * 100
  ) %>% 
  dplyr::select(country_pooled,
                n_plots_total,
                n_lost_26,
                n_lost_45,
                n_lost_85,
                share_lost_26,
                share_lost_45,
                share_lost_85) %>% 
  distinct() %>% 
  mutate(
    country_full = case_when(
      country_pooled == "AT" ~ "Austria",
      country_pooled == "CH" ~ "Switzerland",
      country_pooled == "CZ" ~ "Czech Republic",
      country_pooled == "DE" ~ "Germany",
      country_pooled == "FR" ~ "France",
      country_pooled == "PL" ~ "Poland",
      country_pooled == "SI" ~ "Slovenia",
      country_pooled == "SK" ~ "Slovakia",
      TRUE ~ NA_character_  # Catch anything unexpected
    )
  ) %>%
  arrange(country_full)



#### merge climate suitability tables: richenss and lost plots  ----------------------------
df_out <- full_join(df_stem_suitability_avg_wide,
                    df_suitability_plots) %>% 
  full_join(df_suitability_richness_share)

##### simpler version for MS
df_out_MS <- df_out %>% 
  dplyr::select(country_full,
                #n_plots_total,
                stem_share_rcp26 ,
                stem_share_rcp45 , 
                stem_share_rcp85 ,
                share_lost_26,
                share_lost_45,
                share_lost_85)

##### Full version for Supplement 
df_out_supplement <- df_out %>% 
  dplyr::select(-country_pooled) %>% 
  dplyr::select(country_full,
                richness,  # average number of species per plot/country
                richness26,
                stem_share_rcp26 ,
                richness45,
                stem_share_rcp45 , 
                richness85,
                stem_share_rcp85 ,
                n_plots_total,
                n_lost_26,
                share_lost_26,
                n_lost_45,
                share_lost_45,
                n_lost_85,
                share_lost_85)




### Get Summary row: bottom: Totals : values for richness and lost plots across all values:  -----------
#### Richness & remaining stems 
total_richness_stems_full <- df_suitability_observed_vs_modeled_reg %>% 
  summarize(
    # calculate number of species
    richness = length(unique(acc[current == 1])),
    richness26       = length(unique(acc[current == 1 & rcp26 == 1])),
    richness45       = length(unique(acc[current == 1 & rcp45 == 1])),
    richness85       = length(unique(acc[current == 1 & rcp85 == 1])),
    
    # Total stems currently present
    stems_current = sum(reg_stem_density[current == 1], na.rm = TRUE),
    
    # Total stems that are both currently present AND climatically suitable
    stems_suitable_26 = sum(reg_stem_density[current == 1 & rcp26 == 1], na.rm = TRUE),
    stems_suitable_45 = sum(reg_stem_density[current == 1 & rcp45 == 1], na.rm = TRUE),
    stems_suitable_85 = sum(reg_stem_density[current == 1 & rcp85 == 1], na.rm = TRUE)
    
  ) %>% 
  mutate(country_full = "Total",
         richness26_share = richness26/richness*100,
         richness45_share = richness45/richness*100,
         richness85_share = richness85/richness*100,
         stem_share_rcp26 = stems_suitable_26/stems_current*100,
         stem_share_rcp45 = stems_suitable_45/stems_current*100,
         stem_share_rcp85 = stems_suitable_85/stems_current*100
  ) %>% 
  dplyr::select(country_full,
                richness,
                richness26,
                stem_share_rcp26,
                richness45,
                stem_share_rcp45,
                richness85,
                stem_share_rcp85
  )

# make simpler version for MS
total_stem_dens_MS <- total_richness_stems_full %>% 
  dplyr::select(country_full,
                stem_share_rcp26 ,
                stem_share_rcp45 ,
                stem_share_rcp85 )


#### Summary row: Total loss -----------------------
total_lost_full <- 
  df_suitability_observed_vs_modeled_reg %>%
  group_by(plot) %>%
  mutate(
    any_present = any(current == 1),
    lost26 = if (any_present[1]) all(rcp26[current == 1] == 0) else FALSE,
    lost45 = if (any_present[1]) all(rcp45[current == 1] == 0) else FALSE,
    lost85 = if (any_present[1]) all(rcp85[current == 1] == 0) else FALSE
  ) %>%
  ungroup() %>%
  #group_by(country_pooled) %>%
  summarize(
    n_plots_total = n_distinct(plot),
    n_lost_26 = n_distinct(plot[lost26]),
    n_lost_45 = n_distinct(plot[lost45]),
    n_lost_85 = n_distinct(plot[lost85]),
    share_lost_26 = n_lost_26 / n_plots_total * 100,
    share_lost_45 = n_lost_45 / n_plots_total * 100,
    share_lost_85 = n_lost_85 / n_plots_total * 100
  ) %>% 
  dplyr::select(
    n_plots_total,
    n_lost_26,
    share_lost_26,
    n_lost_45,
    share_lost_45,
    n_lost_85,
    share_lost_85) %>% 
  distinct() 


total_lost_MS <- total_lost_full %>% 
  dplyr::select(
    # n_plots_total,
    share_lost_26,
    share_lost_45,
    share_lost_85) 

##### bind total rows for richness and lost plots: for MS and full for supplement 
# total row for MS:  (only percentages)
# full total rown for supplemnt 
total_row_MS         <- cbind(total_richness_stems_full, total_lost_MS)
total_row_supplement <- cbind(total_richness_stems_full, total_lost_full)


#### Add total rowns to MS and supplement tables: 

# bind the total rows to the original tables: for MS, for Supplement
df_out_MS_with_total         <- bind_rows(df_out_MS, 
                                          total_row_MS) 
df_out_supplement_with_total <- bind_rows(df_out_supplement, 
                                          total_row_supplement)

# format supplement table nicely
df_out_supplement_formatted <- df_out_supplement_with_total %>%
  mutate(
    lost_26 = paste0(n_lost_26, " (", round(share_lost_26, 1), ")"),
    lost_45 = paste0(n_lost_45, " (", round(share_lost_45, 1), ")"),
    lost_85 = paste0(n_lost_85, " (", round(share_lost_85, 1), ")"),
    
    richness = round(richness, 1),
    richness26_fmt = paste0(round(richness26, 1), " (", round(stem_share_rcp26 , 1), ")"),
    richness45_fmt = paste0(round(richness45, 1), " (", round(stem_share_rcp45 , 1), ")"),
    richness85_fmt = paste0(round(richness85, 1), " (", round(stem_share_rcp85 , 1), ")")
  ) %>%
  dplyr::select(
    country_full,
    richness,
    richness26_fmt,
    richness45_fmt,
    richness85_fmt,
    n_plots_total,
    lost_26,
    lost_45,
    lost_85
  )

View(df_out_supplement_formatted)

sjPlot::tab_df(df_out_MS_with_total,
               show.rownames = FALSE,
               file = file.path(public_dir, "tables", "table_MS_clim_suitability_by_country.doc"),
               digits = 1)

sjPlot::tab_df(df_out_supplement_formatted,
               show.rownames = FALSE,
               file = file.path(public_dir, "tables", "table_Supplement_clim_suitability_by_country.doc"),
               digits = 1)


# which species will remain?
# Calculate the presence proportion for each scenario
species_presence_proportion <- df_suitability_observed_vs_modeled %>%
  ungroup() %>%
  dplyr::filter(current == 1) %>%
  dplyr::group_by(acc) %>%
  dplyr::summarise(
    current_count = sum(current == 1, na.rm = T),
    rcp26_count = sum(rcp26 == 1, na.rm = T),
    rcp45_count = sum(rcp45 == 1, na.rm = T),
    rcp85_count = sum(rcp85 == 1, na.rm = T),
    current_proportion = (current_count / 849) * 100,
    rcp26_proportion = (rcp26_count / 849) * 100,
    rcp45_proportion = (rcp45_count / 849) * 100,
    rcp85_proportion = (rcp85_count / 849) * 100
  ) %>%
  arrange(desc(current_count))

# Display the results
View(species_presence_proportion)


# how may plots have ONLY piab??
plots_with_only_piab <- df_suitability_observed_vs_modeled %>%
  dplyr::filter(reg_stem_density > 0) %>%
  group_by(plot) %>%
  summarise(only_piab = all(acc == "piab"), .groups = "drop") %>%
  dplyr::filter(only_piab) %>%
  pull(plot)  

length(plots_with_only_piab)
# 48 from 849, 5.7%

