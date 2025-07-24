
# Process: 
# Get indicators from simulated data
# - read iLand simulated data - simulated by Kilian on 06/11/2024

# Calculate indicators from simulated data
# temporal develpment
#- evaluate they cchange change over time, how thtable they are given the clim cluster

# - inspect the development of 12 iLand landscapes

# simulated data:
# clim_model: 3
# clim_scenario: 4 
# seeds present/not: 2
# repetition: 5
# landscapes: 12




# check data:
# inspect my variables with the simulated ones in year 0/1

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(tidyr)


# Cluster analysis
library(cluster)

# Input data -----------------------------------------------------------
# get simulated data
df_sim <- fread('outTable/df_simulated.csv')
head(df_sim)

# get field data to compare with simulated ones - stem density
df_field      <- fread('outData/veg_density_DBH.csv')
df_indicators <- fread('outData/indicators_for_cluster_analysis.csv') # summarized on plot level

# structural clusters numbers 
df_env_stnd_clust <- fread('rawData/iLand/Cluster_Plots.csv')  # structural and clusters indications for analysis from Kilian

# Plot level vegetation data -----------
# get stem density for juveniles and saplings: to check advanced vs delayed regeneration
# eg is advanced regeneration more stable then delayed ones?

# final tables on site level --------------
df_fin <- fread('outData/indicators_for_cluster_analysis.csv')
df_delayed_advanced <- fread('outTable/df_delayed_advanced.csv') # indication of the delayed, advanced vs other -
# compare the stem desity development across this

# rename clusters from Kilian: env_stnd_clust 
df_sim <- df_sim %>% 
  dplyr::rename(env_stnd_clust  = cluster)

df_delayed_advanced <- df_delayed_advanced %>% 
  mutate(
    adv_delayed = dplyr::recode(adv_delayed, "Other" = "Intermediate"))

# create table to subset the  field data into landscapes
# Create the data frame with the given pairs
df_sites_clusters <- data.frame(
  site = c("23_132", 
           "26_134", "15_133", "17_104", "22_101", "12_151",
           "24_146", "20_116", "12_117", "11_145", "19_160", "25_150"),
  env_stnd_clust  = c("1_1", "1_2", "1_3", "1_4", "2_1", "2_2",
              "2_3", "2_4", "2_5", "3_1", "3_2", "3_3")
)


df_delayed_advanced_sub <- df_delayed_advanced %>% 
  right_join(df_sites_clusters)

# split into delayed vs advanced
df_vegetation_sub <- df_fin %>% 
  ungroup() %>% 
  dplyr::select(site,
               # stem_density,
                sum_stems_juvenile, 
                sum_stems_sapling, 
                sum_stems_mature) %>% 
  mutate(stem_regeneration = sum_stems_juvenile + sum_stems_sapling) %>% 
  mutate(adv_delayed = ifelse(stem_regeneration <= 0, "Delayed", 
                              ifelse(sum_stems_juvenile >= 1000, "Advanced", "Intermediate"))) %>% # Add the 3rd category
  right_join(df_sites_clusters)

length(unique(df_vegetation_sub$site))

# List average field data as input for the landscape level simulation
# "23_132" - "1_1" # 3 vertical classes
# "26_134" - "1_2" # 3 vertical classes
# "15_133" - "1_3"
# "17_104" - "1_4"
# "22_101" - "2_1"
# "12_151" - "2_2"
# "24_146" - "2_3"
# "20_116" - "2_4"
# "12_117" - "2_5"
# "11_145" - "3_1"
# "19_160" - "3_2"
# "25_150" - "3_3"

table(df_vegetation_sub$adv_delayed)


# name the clusters from Kilian: if cluster analysis is run separately, the numbers do not fit! therefore,
# check the naming from plots (boxplot) an rename them to fit
df_sim_class <- df_sim %>%
  mutate(landscape_run = paste(env_stnd_clust, run_nr, sep = '_'),
         seed_comb = paste(seed_flag, seed_sensitivity, sep = "_"))  # yunique run per lanscape scenario



# add delayed vs advanced regeneration indication ------------------------------
# simulated 
df_sim_class <- df_sim_class %>% 
  mutate(unique_sim_run = paste(clim_model, clim_scenario, env_stnd_clust, run_nr, seed_flag, seed_sensitivity, sep = "_")) %>%  # yunique run per lanscape scenario
  right_join(df_delayed_advanced_sub, by = join_by(env_stnd_clust))


## filter field data for selected clusters (simulated landscapes) --------------
df_field_sub <- df_field %>%
  rename(site = cluster) %>%
  right_join(df_sites_clusters, by = join_by(site)) #%>%
  #mutate(clim_cluster = str_sub(landscape, 1, 1),  # add indication of the climatic cluster (1,2,3)
      #   str_cluster = str_sub(landscape, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)


# get structural and climate clusters for all locations:
df_indicators <- df_indicators %>% 
  left_join(df_env_stnd_clust, by = c('site' = 'plot'))



# data check in: check initial values with simulated  ---------------------------------------------------------------------
### select only 12 landscapes form field data --------------------------------------
df_field_ind_sub <- df_indicators %>% 
  dplyr::select(site,  rIVI, richness, stem_density, 
                n_vertical, tmp, prcp, spei3, clay_extract, sand_extract,
                env_stnd_clust, env_cluster, stnd_cluster) %>% 
  right_join(df_sites_clusters) #%>% 

# chcek composite score: try to order variables ----






# Process simulated data --------------------------------------------------------


## Inspect simulated data: ------------------------------------------------------
# # 12 clusters, 
# # 3 climatic models
# # 4 climate scenarios
# # 2 seed added/no
# unique(df_sim2$landscape)  # 12, 3 climatic, 4 strcutural in each climate cluster
# #[1] "1_1" "1_2" "1_3" "1_4" "2_1" "2_2" "2_3" "2_4" "2_5" "3_1" "3_2" "3_3"
# 
# unique(df_sim2$species) # 24
# #[1] "acps" "cabe" "potr" "saca" "quro" "acca" "frex" "ulgl" "lade" "tico" "piab" "fasy" "psme"
# #[14] "soau" "pisy" "bepe" "abal" "alin" "algl" "soar" "coav" "acpl" "rops" "casa"
# 
# 
# unique(df_sim2$clim_model)
# # [1] "ICHEC" "MPI"   "NCC"  
# 
# unique(df_sim2$clim_scenario)
# # [1] "HISTO" "RCP26" "RCP45" "RCP85"
# 
# # external seed: external seed present or not?? TRUE/FALSE
# # years: 1-30
# unique(df_sim2$run_nr)  # 5 repetitions
# 
# 

### Indicators from simulated data : ------------------------------------------------------ 

# richness
# rIVI
# stem density
# vertical classes
# summarize across simulation repetition and across the clim model


#### Species richness ------------------------------------

# Get grouping table : 1380 unique simulations combination
df_simulations_groups <- df_sim_class %>% 
  dplyr::select(year, clim_model, seed_comb,landscape_run, unique_sim_run) %>% 
  distinct() #%>%  # remove duplicated rows
  #mutate(landscape = paste(clim_cluster,str_cluster, sep = "_"))
# nrows with years: 42720
  
  
df_simulations_groups_simple <- df_sim_class %>% 
    dplyr::select(clim_model,seed_comb,landscape_run, unique_sim_run) %>% 
    distinct() #%>%  
  
# Richness
df_richness <- 
  df_sim_class %>% 
  dplyr::filter(count_ha >0) %>% 
  group_by(year,  unique_sim_run) %>%  #clim_modelclim_cluster, 
  summarize(richness = n_distinct(species))



###### Species importance value (relative, from rel_density and rel_BA) ------------------------------------
# get species importance value: from relative density, relative BA
# first calculate the total values per ha, then add it to original table to calculate teh rIVI based on relative dominance
df_sum_cluster <- df_sim_class %>% 
  group_by(year, unique_sim_run) %>%  #clim_modelclim_cluster, 
  summarize(sum_stems = sum(count_ha, na.rm = T),
            sum_BA    = sum(basal_area_m2, na.rm = T)) %>% 
  ungroup()


df_IVI <- 
  df_sim_class %>% 
  # group by spei3es across the levels
  group_by(year, species, unique_sim_run) %>%  #clim_modelclim_cluster, 
  summarize(sp_dens = sum(count_ha),
            sp_BA   = sum(basal_area_m2)) %>%
  #nrow()
  ungroup(.) %>% 
    left_join(df_sum_cluster) %>% # merge sums per whole cluster
  mutate(rel_dens  = sp_dens/sum_stems,
         rel_BA  = sp_BA/sum_BA,
         rIVI = ( rel_dens +rel_BA)/2) %>% # relative species importance value
  #nrow()
  dplyr::select(-sp_dens, -sp_BA,-sum_stems,-sum_BA,-rel_dens, -rel_BA) %>%  # remove unnecessary cols
    
  # filter teh dominant species - highest rIVI
  group_by(year,  unique_sim_run) %>% #clim_model ,
  dplyr::filter(rIVI == max(rIVI)) %>%
  sample_n(1) %>%  # Select a random row if there are ties
  rename(dominant_species      = species) 

  
  
  

##### Structure: Vertical classes & stem density -------------------------------------------------------------
# inspeact vertical classes
df_structure <- df_sim_class %>% 
  dplyr::filter(count_ha >0) %>% 
  group_by(year, unique_sim_run) %>%  #clim_modelclim_cluster, 
  summarize(n_vertical = n_distinct(category),
            stem_density = sum(count_ha)) 


# check lengths
nrow(df_richness)
nrow(df_IVI)  # two rows are missing!
nrow(df_structure)

head(df_richness)
head(df_IVI)
head(df_structure)

length(unique(df_richness$unique_sim_run))
length(unique(df_IVI$unique_sim_run))
length(unique(df_structure$unique_sim_run))


# seems that ICHEC_HISTO_seed_3_2 has not noseed alternative?


# merge df indicators 
df_sim_indicators <- 
  df_simulations_groups %>% # keep all simulations
  left_join(df_richness) %>% 
  left_join(df_IVI) %>%
  left_join(df_structure) %>%
  mutate(env_stnd_clust = substr(landscape_run, 1, 3)) 
  
# landscaepe 3_2 has only 'seed' scenario: no seed would lead to no tree regeneration
# replace NA values by 0 
df_sim_indicators <- df_sim_indicators %>%
  mutate(
    richness = ifelse(is.na(richness), 0, richness),
    stem_density = ifelse(is.na(stem_density), 0, stem_density),
    rIVI = ifelse(is.na(rIVI), 0, rIVI),
    n_vertical = ifelse(is.na(n_vertical), 0, n_vertical)
  )



# add indication of advanced vs delayed vegetation: advanced >1000 juveniles, delayed < 50 stems/ha
df_sim_indicators <- df_sim_indicators %>% 
  left_join(df_delayed_advanced_sub, by = join_by(env_stnd_clust))

# get stem density stats in year 30: keeps seeds in, as it is more realistic scenario that no seeds at all
df_sim_indicators %>% 
 # dplyr::filter(ext_seed == 'seed') %>% 
  dplyr::filter(year == 30) %>% 
  group_by(adv_delayed) %>% 
  summarise(median = median(stem_density),
            iqr = IQR(stem_density),
            mean = mean(stem_density),
            sd = sd(stem_density),
            max = max(stem_density))
  





##  with shaded zones ------------------

# Calculate IQR for each year and adv_delayed group
df_summary <- df_sim_indicators %>%
  dplyr::filter(!str_detect(seed_comb, "noseed")) %>% 
  group_by(year, adv_delayed,seed_comb ) %>%
  summarize(
    median_density = median(stem_density, na.rm = TRUE),
    Q1 = quantile(stem_density, 0.25, na.rm = TRUE),
    Q3 = quantile(stem_density, 0.75, na.rm = TRUE)
  ) %>% 
  mutate(adv_delayed = factor(adv_delayed,
                              levels = c("Delayed",  "Intermediate","Advanced" )))  


df_summary_simpl <- df_sim_indicators %>%
  #dplyr::filter(ext_seed == 'seed') %>%
  dplyr::filter(!str_detect(seed_comb, "noseed")) %>% 
  group_by(year, adv_delayed) %>%
  summarize(
    median_density = median(stem_density, na.rm = TRUE),
    Q1 = quantile(stem_density, 0.25, na.rm = TRUE),
    Q3 = quantile(stem_density, 0.75, na.rm = TRUE)
  )

# Ensure the levels of adv_delayed are in the desired order
df_sim_indicators <- df_sim_indicators %>%
    mutate(adv_delayed = factor(adv_delayed, levels = c("Delayed",  "Intermediate","Advanced" )))

# Plot with median line and IQR ribbon
#p_simulated_stem_dens <- 
 # df_sim_indicators %>% 
colors <- setNames(
  colorRampPalette(c("#A50026", "#FDAE61", "#006837"))(7),
  c("seed_-50", "seed_-25", "seed_-10", "seed_0", "seed_10", "seed_25", "seed_50")
)

# Create the plot
ggplot(df_summary, aes(x = year, y = median_density, 
                       color = seed_comb, 
                       fill = seed_comb)) +
  geom_line(size = 1) +  # Line for the median
  geom_ribbon(aes(ymin = Q1, ymax = Q3), alpha = 0.2, color = NA) + # ribbon for Q1-Q3
  facet_wrap(~ adv_delayed) +  # Facet by adv_delayed 
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  labs(title = "",
       x = "Year",
       y = "Stem Density",
       color = "",
       linetype = "",#Reg. status
       fill = "") +#Reg. status
  theme_classic2() +
  theme(legend.position = 'bottom',
        panel.border = element_rect(color = "black", linewidth = 0.7, fill = NA),
        text = element_text(size = 8),             # Set base text size
        axis.text = element_text(size = 8),        # Axis tick labels
        axis.title = element_text(size = 8),       # Axis titles
        strip.text = element_text(size = 8),       # Facet labels
        legend.text = element_text(size = 8),      # Legend text
        legend.title = element_text(size = 8),     # Legend title
        plot.title = element_text(size = 8)        # Plot title)
)


reg_colors <- c(
  "Delayed" = "#A50026",   # reddish
  "Intermediate"   = "#FDAE61",   # yellowish
  "Advanced"= "#006837"    # green
)

reg_colors_short <- c(
  "Del." = "#A50026",   # reddish
  "Int."   = "#FDAE61",   # yellowish
  "Adv."= "#006837"    # green
)
  

# Manually define lighter versions for 'no seed' scenario
reg_colors_pale <- c(
  "Delayed" = "#F4A6A6",    # pale red
  "Intermediate"   = "#FDD9A0",    # pale yellow-orange
  "Advanced"= "#A6D8A8"     # pale green
)

# Create a combined fill vector based on both ext_seed and adv_delayed
fill_values <- c(
  "seed_Delayed"     = reg_colors[["Delayed"]],
  "seed_Intermediate"       = reg_colors[["Intermediate"]],
  "seed_Advanced"    = reg_colors[["Advanced"]],
  "noseed_Delayed"   = reg_colors_pale[["Delayed"]],
  "noseed_Intermediate"     = reg_colors_pale[["Intermediate"]],
  "noseed_Advanced"  = reg_colors_pale[["Advanced"]]
)


# get simple summary -------------------
df_summary_simpl <- df_summary_simpl %>% 
  # mutate(
  #   adv_delayed = recode(adv_delayed, "Other" = "Intermediate")
  # ) %>% 
  mutate(adv_delayed = factor(adv_delayed, 
                              levels = c("Delayed", "Intermediate", "Advanced")))# %>%  # Ensure correct order
  
p_simulated_stem_dens <- 
   df_sim_indicators %>% 
  dplyr::filter(!str_detect(seed_comb, "noseed")) %>% 
  mutate(adv_delayed = factor(adv_delayed, 
                              levels = c("Delayed", "Intermediate", "Advanced"))) %>%  # Ensure correct order
  ggplot(aes(x = year, y = stem_density/1000, 
               group = adv_delayed , 
               color = adv_delayed)) + 
    # Add IQR ribbon
    geom_ribbon(data = df_summary_simpl, aes(x = year, ymin = Q1/1000, ymax = Q3/1000, 
                                             fill = adv_delayed), 
                alpha = 0.2, inherit.aes = FALSE) +  # Adjust transparency for the ribbon
    # Add median line
    geom_line(data = df_summary_simpl, aes(x = year, y = median_density/1000, 
                                     group = adv_delayed,
                                     color = adv_delayed,
                                     linetype = adv_delayed), 
              linewidth = 1, inherit.aes = FALSE) + 
    facet_grid(.~adv_delayed) +
  scale_color_manual(values = reg_colors) +
  scale_fill_manual(values = reg_colors)+
  labs(title = "",
         x = "Year",
         #y = "Stem Density",
         y = expression("Stem density [1000 n ha"^{-1}*"]"),
         color = "",
         linetype = "",#Reg. status
         fill = "") +#Reg. status
    theme_classic2() +
    theme(legend.position = 'none',
          panel.border = element_rect(color = "black", linewidth = 0.7, fill = NA),
          text = element_text(size = 8),             # Set base text size
          axis.text = element_text(size = 8),        # Axis tick labels
          axis.title = element_text(size = 8),       # Axis titles
          strip.text = element_text(size = 8),       # Facet labels
          legend.text = element_text(size = 8),      # Legend text
          legend.title = element_text(size = 8),     # Legend title
          plot.title = element_text(size = 8)        # Plot title)
    )
  
p_simulated_stem_dens
ggsave(filename = 'outFigs/fig_p_simulated_stem_dens.png', 
       plot = p_simulated_stem_dens, width = 6.5, height = 3, dpi = 300, bg = 'white')


# Sensitivity analysis : range of seed and no seed scenario ---------------------------------------

df_sim_indicators <- df_sim_indicators %>%
  mutate(
    # extract numeric seed level from string like "seed_-25"
    seed_level_num = as.integer(str_extract(seed_comb, "-?\\d+")),
    seed_comb = factor(seed_comb, levels = paste0("seed_", sort(unique(seed_level_num))))
  )

library(data.table)

# Save the table using fwrite()
#fwrite(df_sim_indicators, 
 #      "outData/public/data/df_sim_indicators.csv")


df_sens_plot <- df_sim_indicators %>% 
 # dplyr::filter(!str_detect(seed_comb, "noseed")) %>% 
  dplyr::filter(year %in% 25:30) %>%
  mutate(
    seed_comb = ifelse(is.na(seed_comb), "No seed", as.character(seed_comb)),
    seed_level_num = as.integer(str_extract(seed_comb, "-?\\d+")),
    seed_comb = factor(seed_comb, levels = c("No seed", paste0("seed_", 
                                                               sort(unique(seed_level_num))))),
    adv_delayed = recode_factor(adv_delayed,
                                "Delayed" = "Del.",
                                "Intermediate" = "Int.",
                                "Advanced" = "Adv."),
    adv_delayed = factor(adv_delayed, levels = c("Del.", "Int.", "Adv."))
  )

group_summary <- df_sens_plot %>%
  dplyr::group_by(adv_delayed, seed_comb) %>%
  dplyr::summarize(
    q25 = quantile(stem_density / 1000, 0.25, na.rm = TRUE),
    median_density = median(stem_density / 1000, na.rm = TRUE),
    q75 = quantile(stem_density / 1000, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p_simulated_stem_sensitivity <- group_summary %>%
  ggplot(aes(y = seed_comb, x = median_density,
             xmin = q25, xmax = q75,
             color = adv_delayed)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.3) +
  geom_vline(data = group_summary %>%
               group_by(adv_delayed) %>%
               summarize(median_density = median(median_density), .groups = "drop"),
             aes(xintercept = median_density, color = adv_delayed),
             linetype = "dashed", show.legend = FALSE) +
  scale_color_manual(values = reg_colors_short) +
  theme_classic2() +
  scale_x_continuous(limits = c(0, 13)) +
  scale_y_discrete(labels = ~ str_remove(., "seed_")) +
  labs(y = expression("Change in seeds availability [%]"),
       x = expression("Stem density [" * 1000 *~ n~ha^{-1} * "]"),
       color = "") +
  theme(
    legend.position = 'bottom',
    panel.border = element_rect(color = "black", linewidth = 0.7, fill = NA),
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 8)
  )

p_simulated_stem_sensitivity

ggsave(filename = 'outFigs/p_simulated_sensitivity.png', 
       plot = p_simulated_stem_sensitivity, width = 3, height = 4, dpi = 300, bg = 'white')




