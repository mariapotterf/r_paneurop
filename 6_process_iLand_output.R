
# Process: 
# Get indicators from simulated data
# - read iLand simulated data - simulated by Kilian on 06/11/2024

# Calculate indicators from simulated data
# temporal develpment
#- evaluate they cchange change over time, how thtable they are given the clim cluster


# pre-analysis
# - inspect the development of 12 iLand landscapes

# compare stem density bbetween intial condistion and between clusters


# Process all sites:
# - run k means clustering to classify clim_cluster
# - compare my indicators with simulated data over time
# - clim_clusters: 
# --	Cluster 1 – wet, cold
# --	Cluster 2 – hot, dry
# --	Cluster 3 – in between 1-2

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
df_field     <- fread('outData/veg_density_DBH.csv')
df_indicators <- fread('outData/indicators_for_cluster_analysis.csv') # summarized on plot level

# structural clusters numbers 
df_str_compos_clusters_full <- fread('rawData/iLand/Cluster_Plots.csv')  # structural and clusters indications for analysis from Kilian

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
  mutate(adv_delayed = ifelse(stem_regeneration <= 50, "Delayed", 
                              ifelse(sum_stems_juvenile >= 1000, "Advanced", "Other"))) %>% # Add the 3rd category
  right_join(df_sites_clusters)

length(unique(df_vegetation_sub$site))

# List average field data as input for the landscape level simulation
# "23_132" - "1_1"
# "26_134" - "1_2"
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
  # mutate(clim_class = case_when(
  #   clim_cluster == 1 ~ "wet-warm-clay",  # Cluster 1: wet, cold, clay
  #   clim_cluster == 2 ~ "hot-dry-clay",  # Cluster 2: hot, dry, clay (more sand, less clay, more av.nitro than cluster 3)
  #   clim_cluster == 3 ~ "hot-dry-sand"    # Cluster 3: hot, dry, more sand
  # ))  %>% 
  #mutate(landscape_run = paste(env_stnd_clust, run_nr))  # yunique run per lanscape scenario

str(df_sim_class)
head(df_sim_class)


# get delayed vs advanced regeneration indication ------------------------------
# simulated 
df_sim_class <- df_sim_class %>% 
  mutate(unique_sim_run = paste(clim_model, clim_scenario, ext_seed, env_stnd_clust, run_nr, sep = "_")) %>%  # yunique run per lanscape scenario
  right_join(df_delayed_advanced_sub, by = join_by(env_stnd_clust))


## filter field data for selected clusters (simulated landscapes) --------------
df_field_sub <- df_field %>%
  rename(site = cluster) %>%
  right_join(df_sites_clusters, by = join_by(site)) #%>%
  #mutate(clim_cluster = str_sub(landscape, 1, 1),  # add indication of the climatic cluster (1,2,3)
      #   str_cluster = str_sub(landscape, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)

my_cluster_test = '1_2'
# check per one clutsre! 1_2
simulated_test <- df_sim_class %>% 
  dplyr::filter(env_stnd_clust == my_cluster_test)

df_field_test <- df_field_sub %>% 
  dplyr::filter(env_stnd_clust == my_cluster_test)

View(simulated_test)
View(df_field_test)

# get Landsscape indicators for all of locations:
df_indicators <- df_indicators %>% 
 # rename(site = cluster) %>% 
  left_join(df_sites_clusters) %>%
 # mutate()
  mutate(clim_cluster_spei3 = factor(clim_cluster_spei3)) %>%   # add indication of the strutural cluster (1,2,3,4,5)
  mutate(clim_class = case_when(
    clim_cluster_spei3  == 1 ~ "wet-warm-clay",  # Cluster 1: wet, cold, clay
    clim_cluster_spei3  == 2 ~ "hot-dry-clay",  # Cluster 2: hot, dry, clay (more sand, less clay, more av.nitro than cluster 3)
    clim_cluster_spei3  == 3 ~ "hot-dry-sand"    # Cluster 3: hot, dry, more sand
  ))  %>%
mutate(
  scaled_stem_density = (stem_density - min(stem_density)) / (max(stem_density) - min(stem_density)),
  scaled_n_vertical = (n_vertical - min(n_vertical)) / (max(n_vertical) - min(n_vertical)),
  scaled_richness = (richness - min(richness)) / (max(richness) - min(richness)),
  scaled_rIVI = (rIVI - min(rIVI)) / (max(rIVI) - min(rIVI))
) %>%
  # Calculate composite score: higher values for stem_density, n_vertical, richness; lower values for rIVI
  mutate(
    composite_score = round(scaled_stem_density + scaled_n_vertical + scaled_richness - scaled_rIVI, 2)#,
    # cluster = factor(cluster, levels = cluster[order(composite_score, decreasing = TRUE)])
  ) %>% 
  mutate(composite_class = factor(composite_score)) %>% 
  left_join(df_str_compos_clusters_full, by = c('site' = 'plot'))



# data check in: check initial values with simulated  ---------------------------------------------------------------------
### select only 12 landscapes form field data --------------------------------------
df_field_ind_sub <- df_indicators %>% 
  dplyr::select(site,  rIVI, richness, stem_density, 
                n_vertical, tmp, prcp, spei3, clay_extract, sand_extract,
                #clim_cluster_spei3, clim_class,
                env_stnd_clust, env_cluster, stnd_cluster) %>% 
#  rename(site = cluster) %>% 
  right_join(df_sites_clusters) %>% 
  mutate(clim_cluster_spei3 = factor(clim_cluster_spei3)) #%>%  # add indication of the strutural cluster (1,2,3,4,5)
 # rename(landscape = cluster)

# chcek composite score: try to order variables ----

# Standardize the variables
df_field_ind_sub <- df_field_ind_sub %>%
  mutate(
    scaled_stem_density = (stem_density - min(stem_density)) / (max(stem_density) - min(stem_density)),
    scaled_n_vertical = (n_vertical - min(n_vertical)) / (max(n_vertical) - min(n_vertical)),
    scaled_richness = (richness - min(richness)) / (max(richness) - min(richness)),
    scaled_rIVI = (rIVI - min(rIVI)) / (max(rIVI) - min(rIVI))
  ) %>%
  # Calculate composite score: higher values for stem_density, n_vertical, richness; lower values for rIVI
  mutate(
    composite_score = round(scaled_stem_density + scaled_n_vertical + scaled_richness - scaled_rIVI, 2)#,
   # cluster = factor(cluster, levels = cluster[order(composite_score, decreasing = TRUE)])
  ) %>% 
  mutate(composite_class = factor(composite_score)) %>% 
  mutate(clim_class = case_when(
    clim_cluster_spei3  == 1 ~ "wet-warm-clay",  # Cluster 1: wet, cold, clay
    clim_cluster_spei3  == 2 ~ "hot-dry-clay",  # Cluster 2: hot, dry, clay (more sand, less clay, more av.nitro than cluster 3)
    clim_cluster_spei3  == 3 ~ "hot-dry-sand"    # Cluster 3: hot, dry, more sand
  ))  #%>% 





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
  dplyr::select(year, clim_model, ext_seed,landscape_run, unique_sim_run) %>% 
  distinct() #%>%  # remove duplicated rows
  #mutate(landscape = paste(clim_cluster,str_cluster, sep = "_"))
# nrows with years: 42720
  
  
df_simulations_groups_simple <- df_sim_class %>% 
    dplyr::select(clim_model,ext_seed,landscape_run, unique_sim_run) %>% 
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
  mutate(landscape = substr(landscape_run, 1, 3))
  
# landscaepe 3_2 has only 'seed' scenario: no seed would lead to no tree regeneration
# replace NA values by 0 
df_sim_indicators <- df_sim_indicators %>%
  mutate(
    richness = ifelse(is.na(richness), 0, richness),
    stem_density = ifelse(is.na(stem_density), 0, stem_density),
    rIVI = ifelse(is.na(rIVI), 0, rIVI),
    n_vertical = ifelse(is.na(n_vertical), 0, n_vertical)
  )

# 
# get a composite index to order the landscapes
df_sim_indicators <- df_sim_indicators %>% 
  mutate(
    scaled_stem_density = (stem_density - min(stem_density)) / (max(stem_density) - min(stem_density)),
    scaled_n_vertical = (n_vertical - min(n_vertical)) / (max(n_vertical) - min(n_vertical)),
    scaled_richness = (richness - min(richness)) / (max(richness) - min(richness)),
    scaled_rIVI = (rIVI - min(rIVI)) / (max(rIVI) - min(rIVI))
  ) %>%
  # Calculate composite score: higher values for stem_density, n_vertical, richness; lower values for rIVI
  mutate(
    composite_score = round(scaled_stem_density + scaled_n_vertical + scaled_richness - scaled_rIVI, 2)) 


# add indication of advanced vs delayed vegetation: advanced >1000 juveniles, delayed < 50 stems/ha
df_sim_indicators <- df_sim_indicators %>% 
  left_join(df_delayed_advanced_sub, by = join_by(landscape))

# get stem density stats in year 30
df_sim_indicators %>% 
  dplyr::filter(ext_seed == 'seed') %>% 
  dplyr::filter(year == 30) %>% 
  group_by(adv_delayed) %>% 
  summarise(median = median(stem_density),
            iqr = IQR(stem_density),
            mean = mean(stem_density),
            sd = sd(stem_density),
            max = max(stem_density))
  


# Plot: differentiate only by the advanced vs delayed ---------------------------
# only stem density
df_sim_indicators %>% 
  ggplot(aes(x = year,
             y = stem_density,
             group= unique_sim_run,
             color =adv_delayed )) + 
  geom_line(alpha = 0.1) +
  facet_grid(. ~adv_delayed) +
  stat_summary(
    aes(group = adv_delayed,
        color = adv_delayed),          # Group by `adv_delayed` for the median line
    fun = median,                      # Calculate the median
    geom = "line",                     # Draw the median line
   # color = "black",                   # Set the color of the median line (e.g., black)
    linewidth = 1,                          # Set the line thickness
    linetype = "dashed"                # Optional: make the line dashed for emphasis
  ) +
  theme_classic()



##  with shaded zones ------------------

# Calculate IQR for each year and adv_delayed group
df_summary <- df_sim_indicators %>%
  #dplyr::filter(ext_seed == 'seed') %>% 
  group_by(year, adv_delayed,ext_seed) %>%
  summarize(
    median_density = median(stem_density, na.rm = TRUE),
    Q1 = quantile(stem_density, 0.25, na.rm = TRUE),
    Q3 = quantile(stem_density, 0.75, na.rm = TRUE)
  )

df_summary_simpl <- df_sim_indicators %>%
  dplyr::filter(ext_seed == 'seed') %>% 
  group_by(year, adv_delayed) %>%
  summarize(
    median_density = median(stem_density, na.rm = TRUE),
    Q1 = quantile(stem_density, 0.25, na.rm = TRUE),
    Q3 = quantile(stem_density, 0.75, na.rm = TRUE)
  )



df_summary <- df_summary %>% 
  mutate(adv_delayed = factor(adv_delayed,
                              levels = c("Delayed",  "Other","Advanced" )))  #%>% 
  
# Ensure the levels of adv_delayed are in the desired order
df_sim_indicators <- df_sim_indicators %>%
  
  mutate(adv_delayed = factor(adv_delayed, levels = c("Delayed",  "Other","Advanced" )))

# Plot with median line and IQR ribbon
#p_simulated_stem_dens <- 
 # df_sim_indicators %>% 

# Create the plot
ggplot(df_summary, aes(x = year, y = median_density, color = ext_seed, fill = ext_seed)) +
  geom_line(size = 1) +  # Line for the median
  geom_ribbon(aes(ymin = Q1, ymax = Q3), alpha = 0.2) +  # Ribbon for Q1-Q3
  facet_wrap(~ adv_delayed) +  # Facet by adv_delayed 
 # scale_color_manual(values = c("#A50026", 
 #                               "#FDAE61",
 #                               "#006837")) +
 # scale_fill_manual(values = c("#A50026", 
 #                               "#FDAE61",
 #                               "#006837")) +
 
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

  
  p_simulated_stem_dens <- 
   df_sim_indicators %>% 
    dplyr::filter(ext_seed == 'seed') %>%  # seelct only seeds scenario - more realistic than no seed
    ggplot(aes(x = year, y = stem_density, 
               group = adv_delayed , 
               color = adv_delayed)) + 
    # Add IQR ribbon
    geom_ribbon(data = df_summary_simpl, aes(x = year, ymin = Q1, ymax = Q3, fill = adv_delayed), 
                alpha = 0.2, inherit.aes = FALSE) +  # Adjust transparency for the ribbon
    # Add median line
    geom_line(data = df_summary_simpl, aes(x = year, y = median_density, 
                                     group = adv_delayed,
                                     color = adv_delayed,
                                     linetype = adv_delayed), 
              linewidth = 1, inherit.aes = FALSE) + 
    facet_grid(.~adv_delayed) +
  scale_color_manual(values = c("#A50026", 
                                "#FDAE61",
                                "#006837")) +
    scale_fill_manual(values = c("#A50026", 
                                 "#FDAE61",
                                 "#006837")) +
    
    labs(title = "",
         x = "Year",
         y = "Stem Density",
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

## Evaluate initial state : for all clusters   -------------------------
# filter initial state: year == 0
df_sim_indicators0 <- df_sim_indicators %>% 
  ungroup() %>% 
  dplyr::filter(year == 0 ) %>% #& clim_scenario == "HISTO" %>% 
  dplyr::select(landscape,  rIVI, richness, stem_density, n_vertical, unique_sim_run, ext_seed ) #%>% 
  
# variation at the end of simulation run
df_sim_indicators_end <- df_sim_indicators %>% 
  ungroup() %>% 
  dplyr::filter(year %in% 25:30 ) %>% #& clim_scenario == "HISTO" %>% 
  dplyr::select(landscape,  rIVI, richness, stem_density, n_vertical, unique_sim_run, ext_seed ) #%>% 


# Calculate median and IQR to merge simulated and field data -------------------

# show stem density per str and clim clusters 
# Calculate median and IQR for stem_density grouped by env_stnd_clust
field_stem_density_summary <- df_indicators %>%
  group_by(env_stnd_clust) %>%
  summarise(
    median_stem_density = median(stem_density, na.rm = TRUE),
    iqr_lower = quantile(stem_density, 0.25, na.rm = TRUE),
    iqr_upper = quantile(stem_density, 0.75, na.rm = TRUE)
  )

# Create the bar plot with error bars
ggplot(field_stem_density_summary, aes(x = env_stnd_clust, y = median_stem_density)) +
  geom_bar(stat = "identity", fill = "grey", color = "grey") +
  geom_errorbar(aes(ymin = iqr_lower, ymax = iqr_upper), width = 0.2) +
  labs(
    x = "Environmental Cluster",
    y = "Stem Density (median ± IQR)",
    title = "Median Stem Density by Environmental Cluster"
  ) +
  theme_minimal()


# observed: ---------------------------------------------------------------
simulated_stem_density_summary <- df_sim_indicators %>%
  dplyr::filter(ext_seed  == 'noseed') %>% 
  dplyr::filter(year %in% c(1)) %>% 
  rename( env_stnd_clust = landscape ) %>% 
  group_by(env_stnd_clust) %>%
  summarise(
    median_stem_density = median(stem_density, na.rm = TRUE),
    iqr_lower = quantile(stem_density, 0.25, na.rm = TRUE),
    iqr_upper = quantile(stem_density, 0.75, na.rm = TRUE)
  )

# Create the bar plot with error bars
ggplot(simulated_stem_density_summary, aes(x = env_stnd_clust, y = median_stem_density)) +
  geom_bar(stat = "identity", fill = "grey", color = "grey") +
  geom_errorbar(aes(ymin = iqr_lower, ymax = iqr_upper), width = 0.2) +
  labs(
    x = "Environmental Cluster",
    y = "Stem Density (median ± IQR)",
    title = "Median Stem Density by Environmental Cluster"
  ) +
  theme_minimal()





###### merge field data with simulated data in year 0 
df_compare0 <-  # env_stnd_clust 
  df_indicators %>% 
  #df_field_ind_sub %>% 
  left_join(df_sim_indicators0, by = c("env_stnd_clust" = "landscape" ), suffix = c("_field", "_simul")) %>%  #by = c( "cluster" = "landscape"), 
  dplyr::select(landscape, site, ext_seed, ends_with("_field"), ends_with("_simul")) %>% 
 na.omit() #%>% # remove empty one


table(df_compare0$richness_field)
table(df_compare0$richness_simul)

plot(df_compare0$richness_field, df_compare0$richness_simul)



# Check the updated data with the ordered factor levels
df_field_ind_sub %>% 
  arrange(desc(composite_score)) %>% 
  dplyr::select(cluster, stem_density, n_vertical, richness, rIVI, composite_score)

plot(df_field_ind_sub$composite_score, 
     df_field_ind_sub$rIVI)
# PLOT simple:  fill in all data: field and simulated, add time to field ata yyear = 1 ---------------
# subset indicators for field data
# ffor simulated data - time range
# plot over one plot

df_indicators_field <- df_indicators %>% 
  dplyr::select(site,  rIVI, richness, stem_density, n_vertical) %>%
  mutate(year = 1) %>% 
  mutate(type = 'field',
         landscape_run = 'field') 


df_simulated <- df_sim_indicators %>% 
  ungroup() %>%  
  mutate(type = 'simulated') %>% 
  dplyr::select(landscape,  rIVI, richness, stem_density, n_vertical, year,   type, landscape_run)
 

head(df_indicators_field)

head(df_simulated)


df_merge <- rbind(df_simulated, df_indicators_field)


df_merge %>% 
  ggplot(aes(x = year,
             y = stem_density,
             group = landscape_run)) + 
  geom_line()



# Summary tables -----------------------------------------------------------------

# Step 1: Calculate median and IQR for field data (year = 1)
field_summary <- df_merge %>%
  dplyr::filter(type == 'field') %>% 
  group_by(clim_class) %>% 
 # filter(year == 1) %>%
  summarise(
  
    median_density = median(stem_density),
    IQR_lower = quantile(stem_density, 0.25),
    IQR_upper = quantile(stem_density, 0.75)
  ) %>% 
  mutate(  year = 1)

(field_summary)

# Step 2: Calculate median and IQR for simulated data (year = 30)
sim_summary <- df_merge %>%
  dplyr::filter(type == 'simulated') %>% 
  group_by(clim_class) %>% 
  filter(year == 25) %>%
  summarise(
    year = 25,  # Add the year column for plotting
    median_density = median(stem_density),
    IQR_lower = quantile(stem_density, 0.25),
    IQR_upper = quantile(stem_density, 0.75)
  )

(sim_summary)

df_merge %>% 
  dplyr::filter(type == 'simulated') %>% 
  ggplot() +
  geom_line(aes(x = year, y = stem_density, group = landscape_run), alpha = 0.2) +
  
  # Error bar for field data (year = 1)
  geom_errorbar(data = field_summary, aes(
    x = year,       
    ymin = IQR_lower,
    ymax = IQR_upper,
    color = clim_class
  ), width = 0.2, position = position_dodge(width = 0.5)) +  # Reduce the dodge width
  
  geom_point(data = field_summary, aes(x = year, y = median_density, color = clim_class), 
             size = 3, position = position_dodge(width = 0.5)) +  # Consistent dodge for points
  
  # Error bar for simulated data (year = 30)
  geom_errorbar(data = sim_summary, aes(
    x = year,       
    ymin = IQR_lower,
    ymax = IQR_upper,
    color = clim_class
  ), width = 0.2, position = position_dodge(width = 0.5)) +  # Same dodge width
  
  geom_point(data = sim_summary, aes(x = year, y = median_density, color = clim_class), 
             size = 3, position = position_dodge(width = 0.5)) +  # Consistent dodge for points
  
  scale_color_manual(values = c(
    "wet-warm-clay" = "blue",         
    "hot-dry-sand" = "red",           
    "hot-dry-clay" = "orange"         
  )) +
  
  labs(x = "Year", y = "Stem Density") +
  theme_classic2()


# get structural cluster analysis -------------------------------------------------
# for data 25:30 years after simulation\
head(df_sim_indicators)


# Simulated: run structural cluster anaysis, for each clim cluster separately --------------

# Subset the data
# "dominant_species", 
# note that disturbace characteristica are missing from here!

# select last decade
data_subset_str <- df_sim_indicators %>% 
  # select only one scenario (climate is missing but that is ok)
  dplyr::filter(clim_model  == 'ICHEC' & ext_seed == 'noseed') %>% 
  dplyr::filter(year %in% 25:30)

# keep original data, to have indicators of original cluster numbers
df_origin_values_sankey <- data_subset_str

data_subset_str  <- data_subset_str[, c("dominant_species", 
                                          "rIVI", "richness",  
                               "stem_density", "n_vertical", 
                               "clim_cluster")]



# Convert dominant_species to factor
data_subset_str$dominant_species <- as.factor(data_subset_str$dominant_species)

# Normalize the quantitative variables for clustering
data_normalized <- data_subset_str %>%
  mutate(across(c(rIVI, richness, stem_density, n_vertical), scale)) %>% 
  as.data.frame()

# Check for missing species and their indices
empty_species_indices <- which(is.na(data_subset_str$dominant_species) | data_subset_str$dominant_species == "")


# Dummy code the categorical variable
data_dummied <- model.matrix(~dominant_species - 1, data = data_normalized)
data_dummied <- as.data.frame(data_dummied)

# Create a matrix of NA with the same number of columns as data_dummied
if (length(empty_species_indices) > 0) {
  na_rows <- matrix(NA, nrow = length(empty_species_indices), ncol = ncol(data_dummied))
  colnames(na_rows) <- colnames(data_dummied)
  
  # Combine the dummy coded categorical variable with NA rows
  data_dummied <- rbind(data_dummied, na_rows)
}

# Combine the dummy coded categorical variable with the normalized quantitative variables
data_for_clustering <- cbind(data_dummied, data_normalized[, -1])

# avoid dummy variables, hard to read
# split table in 3 clusters to run separately structurral cluster analysis:
df1 <- data_for_clustering[, -1] %>%  dplyr::filter(clim_cluster == 1)
df2 <- data_for_clustering[, -1] %>%  dplyr::filter(clim_cluster == 2)
df3 <- data_for_clustering[, -1] %>%  dplyr::filter(clim_cluster == 3)

#
# split original table in 3 clusters to run separately structurral cluster analysis:
df1_full <- df_origin_values_sankey %>%  dplyr::filter(clim_cluster == 1)
df2_full <- df_origin_values_sankey %>%  dplyr::filter(clim_cluster == 2)
df3_full <- df_origin_values_sankey %>%  dplyr::filter(clim_cluster == 3)


# remove last columns so they are not part of teh PCA aalysis
df1 <- df1 %>% dplyr::select(-clim_cluster)
df2 <- df2 %>% dplyr::select(-clim_cluster)
df3 <- df3 %>% dplyr::select(-clim_cluster)

nrow(df1)
nrow(df2)
nrow(df3)

# Find which number of clusters is teh best
# Perform K-means clustering for different values of k

# Function to perform K-means clustering and return silhouette widths
find_best_k <- function(data, max_clusters = 10) {
  sil_width <- numeric(max_clusters)
  
  for (k in 2:max_clusters) {
    kmeans_result <- kmeans(data, centers = k, nstart = 25)
    sil <- silhouette(kmeans_result$cluster, dist(data))
    sil_width[k] <- mean(sil[, 3])
  }
  
  # Return the silhouette widths
  return(sil_width)
}

# Find the best number of clusters for each data frame
sil_width_df1 <- find_best_k(df1)
sil_width_df2 <- find_best_k(df2)
sil_width_df3 <- find_best_k(df3)

sil_width_all = cbind(sil_width_df1,
                      sil_width_df2,
                      sil_width_df3)
sil_width_all <- as.data.frame(sil_width_all)

fwrite(sil_width_all, 'outTable/sil_width_test1.csv')
#fwrite(df_compare0,       'outTable/compare_field_sim_12_lands_start.csv')

# Plot Silhouette widths for each data frame
par(mfrow = c(3, 1))  # Arrange plots in a single column

# Plot for df1
plot(1:length(sil_width_all$sil_width_df1), sil_width_df1, type = "b", 
     xlab = "Number of clusters", ylab = "Average Silhouette width", 
     main = "Silhouette Analysis for df1", ylim = c(0, max(c(sil_width_df1, sil_width_df2, sil_width_df3), na.rm = TRUE)))

# Plot for df2
plot(1:length(sil_width_df2), sil_width_df2, type = "b", 
     xlab = "Number of clusters", ylab = "Average Silhouette width", 
     main = "Silhouette Analysis for df2", ylim = c(0, max(c(sil_width_df1, sil_width_df2, sil_width_df3), na.rm = TRUE)))

# Plot for df3
plot(1:length(sil_width_df3), sil_width_df3, type = "b", 
     xlab = "Number of clusters", ylab = "Average Silhouette width", 
     main = "Silhouette Analysis for df3", ylim = c(0, max(c(sil_width_df1, sil_width_df2, sil_width_df3), na.rm = TRUE)))


best_k_df1 <- 2 #which.max(sil_width_df1) # 2
best_k_df2 <- 2 # 
#which.max(sil_width_df2) # 2
best_k_df3 <- 2 # which.max(sil_width_df3) # 3

best_k_df1 <- which.max(sil_width_df1) # 2
best_k_df2 <- which.max(sil_width_df2) # 2
best_k_df3 <- which.max(sil_width_df3) # 3
(best_k_df1)
(best_k_df2)
(best_k_df3)

# Determine the optimal number of clusters
#optimal_k_clim <- best_k_df1  # from Kilian's study

# Set the optimal number of clusters
optimal_k_clim <- 2  # from visual setup

# List of data frames
data_frames <- list(df1 = df1, df2 = df2, df3 = df3)
#data_frames_full <- list(df1_full = df1_full, df2_full = df2_full, df3_full = df3_full)


#set.seed(3)
#kmeans_result <- kmeans(df3, centers = optimal_k_clim, nstart = 25)

# Loop through each data frame, perform K-means, and store the results
# Loop through each data frame, perform K-means, and store the results
for (i in seq_along(data_frames)) {
  df_name <- names(data_frames)[i]  # Get the name of the current data frame
  (print(df_name))
  print(i)
  # Get the current data frame
  current_df <- data_frames[[df_name]]
  
  # Perform K-means clustering with the optimal number of clusters
  set.seed(3)
  kmeans_result <- kmeans(current_df, centers = optimal_k_clim, nstart = 25)
  
  # Create the str_cluster format as "run_index_cluster"
  print(paste(i, kmeans_result$cluster, sep = "_"))
  data_frames[[df_name]]$str_cluster_end <- paste(i, kmeans_result$cluster, sep = "_")
  
}

data_frames[["df1"]]$clim_cluster <- 1
data_frames[["df2"]]$clim_cluster <- 2
data_frames[["df3"]]$clim_cluster <- 3

df1_full$str_cluster_end <- data_frames[["df1"]]$str_cluster_end
df2_full$str_cluster_end <- data_frames[["df2"]]$str_cluster_end
df3_full$str_cluster_end <- data_frames[["df3"]]$str_cluster_end



# Add cluster indication to original data
df_str_clustering_out <- bind_rows(df1_full,df2_full,df3_full)

table(df_str_clustering_out$clim_cluster, df_str_clustering_out$str_cluster_end )



# make sankey plot -----------------------------------------------
library(networkD3)

# Prepare data for the Sankey plot
sankey_data <- df_str_clustering_out %>%
  group_by(clim_cluster, str_cluster, str_cluster_end) %>%
  summarise(count = n(), .groups = 'drop')  # Count occurrences of each cluster pair

# Create unique nodes
nodes <- data.frame(name = unique(c(sankey_data$str_cluster, sankey_data$str_cluster_end)))

# Create links data frame
links <- sankey_data %>%
  mutate(source = match(str_cluster, nodes$name) - 1,  # Match the str_cluster to node indices
         target = match(str_cluster_end, nodes$name) - 1) %>%  # Match the str_cluster_end to node indices
  select(source, target, value = count)  # Select the necessary columns

# Create the Sankey plot
sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target", Value = "value",
              NodeID = "name", units = "TWh", fontSize = 12, nodeWidth = 30)


# for individual groups

# Prepare a list to store plots
sankey_plots <- list()

# Loop through each unique clim_cluster
for (cluster in unique(df_str_clustering_out$clim_cluster)) {
  
  # Prepare data for the Sankey plot for the current clim_cluster
  sankey_data <- df_str_clustering_out %>%
    filter(clim_cluster == cluster) %>%
    group_by(str_cluster, str_cluster_end) %>%
    summarise(count = n(), .groups = 'drop')  # Count occurrences of each cluster pair
  
  # Create unique nodes
  nodes <- data.frame(name = unique(c(sankey_data$str_cluster, sankey_data$str_cluster_end)))
  
  # Create links data frame
  links <- sankey_data %>%
    mutate(source = match(str_cluster, nodes$name) - 1,  # Match the str_cluster to node indices
           target = match(str_cluster_end, nodes$name) - 1) %>%  # Match the str_cluster_end to node indices
    select(source, target, value = count)  # Select the necessary columns
  
  # Create the Sankey plot for the current clim_cluster
  sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target", 
                               Value = "value", NodeID = "name", units = "TWh", fontSize = 12, nodeWidth = 30)
  
  # Store the plot in the list
  sankey_plots[[as.character(cluster)]] <- sankey_plot
}

# To view a specific plot, e.g., for the first clim_cluster
sankey_plots[[2]]  # Adjust the index as needed




# Make a function: 
# Function to create a plot for a given y-variable
plot_variable <- function(df_merge, y_var) {
  
  # Step 1: Calculate median and IQR for field data (year = 1)
  field_summary <- df_merge %>%
    dplyr::filter(type == 'field') %>%
    group_by(clim_class) %>%
    summarise(
      median_value = median(!!sym(y_var)),
      IQR_lower = quantile(!!sym(y_var), 0.25),
      IQR_upper = quantile(!!sym(y_var), 0.75)
    ) %>%
    mutate(year = 1)
  
  # Step 2: Calculate median and IQR for simulated data (year = 25)
  sim_summary <- df_merge %>%
    dplyr::filter(type == 'simulated') %>%
    group_by(clim_class) %>%
    filter(year == 20:30) %>%
    summarise(
      median_value = median(!!sym(y_var)),
      IQR_lower = quantile(!!sym(y_var), 0.25),
      IQR_upper = quantile(!!sym(y_var), 0.75)
    ) %>%
    mutate(year = 25)
  
  # Step 3: Create the plot
  plot <- df_merge %>%
    dplyr::filter(type == 'simulated') %>%
    ggplot() +
    geom_line(aes(x = year, y = !!sym(y_var), group = landscape_run), alpha = 0.2) +
    
    # Error bars for field data (year = 1)
    geom_errorbar(data = field_summary, aes(
      x = year,       
      ymin = IQR_lower,
      ymax = IQR_upper,
      color = clim_class
    ), width = 0.2, position = position_dodge(width = 0.5)) +
    
    geom_point(data = field_summary, aes(x = year, y = median_value, color = clim_class), 
               size = 3, position = position_dodge(width = 0.5)) +
    
    # Error bars for simulated data (year = 25)
    geom_errorbar(data = sim_summary, aes(
      x = year,       
      ymin = IQR_lower,
      ymax = IQR_upper,
      color = clim_class
    ), width = 0.2, position = position_dodge(width = 0.5)) +
    
    geom_point(data = sim_summary, aes(x = year, y = median_value, color = clim_class), 
               size = 3, position = position_dodge(width = 0.5)) +
    
    scale_color_manual(values = c(
      "wet-warm-clay" = "blue",         
      "hot-dry-sand" = "red",           
      "hot-dry-clay" = "orange"         
    )) +
    
    labs(x = "Year", y = y_var) +  # Dynamic y label based on the variable name
    theme_classic2()
  
  return(plot)
}


# Define the variables you want to loop over
y_variables <- c("stem_density", "rIVI", "richness", "n_vertical")

# Loop over each variable and generate the plots
plots <- list()  # Store plots in a list
for (y_var in y_variables) {
  plot <- plot_variable(df_merge, y_var)
  print(plot)  # Display each plot
  plots[[y_var]] <- plot  # Store plot in a list
}

# Print all plots
p1 <- plots[[1]]
p2 <- plots[[2]]
p3 <- plots[[3]]
p4 <- plots[[4]]
#p5 <- plots[[5]]

windows()
ggarrange(p1, p2, p3, p4,common.legend = TRUE , ncol = 2, nrow = 2)


# Add an overarching title and subtitle to the ggarrange layout
# Use ggarrange to combine the plots
combined_plot <- ggarrange(p1, p2, p3, p4, 
                           ncol = 2, nrow = 2, 
                           common.legend = TRUE, legend = "top") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),  # Title format
    plot.subtitle = element_text(size = 10, hjust = 0.5),  # Subtitle format
    axis.title.y = element_text(size = 8),  # Y-axis labels
    axis.text.y = element_text(size = 8),    # Y-axis text
    axis.text.x = element_text(size = 8),    # X-axis text
    legend.title = element_text(size = 10),  # Legend title
    legend.text = element_text(size = 8)     # Legend text size
  )

# Add title and subtitle using annotate_figure
p_field_vs_simulated <- annotate_figure(combined_plot,
                top = text_grob("Development of Indices Over Time by Climatic-Environmental Clusters", 
                                color = "black", face = "bold", size = 12),
                bottom = text_grob("Error bars represent median and IQR. Field data: year 1, simulated data: year 25 (median-IQR for years 20-30)", 
                                   color = "black", size = 10)
)

ggsave(filename = 'outFigs/fig_field_vs_simulated.png', 
       plot = p_field_vs_simulated, width = 7, height = 7, dpi = 300, bg = 'white')


# Plot: indicators development oevr time -------------------------------------

# df_wide <- df_sim_indicators %>%
#   dplyr::select(year, cluster, clim_scenario, ext_seed, richness, dominant_species, rIVI, n_vertical, stem_density, clim_cluster,str_cluster) %>% 
#   pivot_wider(
#     names_from = !clim_scenario,clim_cluster, str_cluster),
#     values_from = c(richness, dominant_species, rIVI, n_vertical, stem_density)
#   )

p_noseed <-df_sim_indicators  %>% 
  dplyr::filter(ext_seed == 'noseed') %>% 
  ggplot(aes(x = year,
             y = rIVI,
             group = cluster,
             color = str_cluster)) +
  geom_line(alpha = .5) + 
  #geom_point() +
  facet_grid(. ~clim_cluster)

p_seed <-df_sim_indicators  %>% 
  dplyr::filter(ext_seed == 'seed') %>% 
  ggplot(aes(x = year,
             y = rIVI,
             group = cluster,
             color = str_cluster)) +
  geom_line(alpha = .5) + 
  #geom_point() +
  facet_grid(. ~clim_cluster)

ggarrange(p_noseed, p_seed, common.legend = TRUE, labels = c('no seed', 'seeds'))



# PLot --------------------------------------------------------------------------
# compare all field indicators with simulated indicators in time 0 and year 30
# make a common table

# Calculate the average values of numeric columns across ext_seed scenarios
df_sim_indicators_avg <- df_sim_indicators %>%
  group_by(year, cluster, clim_scenario, dominant_species, clim_cluster , str_cluster) %>%  # Group by the relevant columns excluding ext_seed
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%  # Calculate the mean for numeric columns
  ungroup()  # Remove grouping


df_sim_indicators_avg_sub <- df_sim_indicators %>%
  group_by(year, cluster, clim_scenario, dominant_species, clim_cluster , str_cluster) %>%  # Group by the relevant columns excluding ext_seed
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')  %>%  
  # Calculate the mean for numeric columns
  ungroup(.) %>% 
  dplyr::select(year, cluster, clim_cluster,clim_scenario , dominant_species, rIVI, richness , stem_density, n_vertical) %>% 
  rename(site = cluster)  %>% 
  mutate(source = 'simulated')

# Create plot overtime, all scenarios, ad on top the range of teh fiel data (corresponds to year  0/1)
df_sim_indicators_avg %>% 
  dplyr::filter(clim_scenario == 'HISTO') %>% 
 # View()
  ggplot(aes(x = year,y = richness, color = site)) +
  geom_jitter() +
  facet_grid(.~clim_cluster)



##### make a simple plot: just compare the ranges -------------------------------
# need to merge simulated and field data - complete teh necessary columns!
df_field_ind_sub <- 
  df_indicators %>% 
  dplyr::select(site, clim_cluster_test, dominant_species,rIVI, richness, stem_density, n_vertical) %>% 
  rename(clim_cluster = clim_cluster_test) %>% 
  mutate(year = 0,
         clim_scenario  = "HISTO") %>% 
  dplyr::select(year, site, clim_cluster, clim_scenario, dominant_species, rIVI, richness, stem_density, n_vertical) %>% 
  mutate(source = 'observed')

 
# merge tables into one
df_indi_merged <- rbind(df_field_ind_sub,df_sim_indicators_avg_sub)

df_indi_merged <- df_indi_merged %>% 
  mutate(clim_class = case_when(clim_cluster  == '1' ~ 'wet-cold',
                                clim_cluster  == '2'  ~ 'hot-dry',
                                clim_cluster  == '3'  ~ 'medium',
                                TRUE ~ NA_character_
                                )) %>% 
  mutate(clim_class = factor(clim_class, levels = c('wet-cold','medium', 'hot-dry' )),
         clim_cluster = factor(clim_cluster, levels = c('1','3', '2' )))

# reorder and rename clim cluster: 
# 1 - cold, wet
# 2 - hot, dry
# 3 - medium, in between 1-2


# Define custom color palette mimicking RdYlBu for discrete values
custom_palette <- c("wet-cold" = "#4575b4", "medium" = "#ffffbf", "hot-dry" = "#d73027")

# 
# Create the ggplot: -------------------------
# compare the fiel ata with medians over 30 years! (or select just the year ==30??)


# simplified:
# Function to create the plot
create_plot <- function(data, y_var, y_label) {
  ggplot(data, aes(x = source, y = !!sym(y_var), fill = clim_class)) +
    geom_violin(position = position_dodge(width = 0.5)) +
    
    stat_summary(
      fun = median,
      fun.min = function(x) { quantile(x, probs = 0.25) },
      fun.max = function(x) { quantile(x, probs = 0.75) },
      geom = "pointrange",
      position = position_dodge(width = 0.5),
      size = 1  # Adjust point and line size
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      color = "white",
      size = 3.5,  # Adjust point size
      position = position_dodge(width = 0.5)
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      aes(color = clim_class),
      size = 2.5,  # Adjust point size
      position = position_dodge(width = 0.5)
    ) +
    scale_fill_manual(values = custom_palette) +  # Use the custom discrete palette
    scale_color_manual(values = custom_palette) +  # Use the custom discrete palette
    labs(
      x = "Source",
      y = y_label,
      color = "Climate Cluster",
      fill = "Climate Cluster"
    ) +
    theme_bw()
}

# Filter the data once
filtered_data <- df_indi_merged %>% 
  dplyr::filter(clim_scenario == "HISTO")

# Create plots
p_IVI <- create_plot(filtered_data, "rIVI", "rIVI")
p_rich <- create_plot(filtered_data, "richness", "Richness")
p_stem_dens <- create_plot(filtered_data, "stem_density", "Stem Density")
p_vert <- create_plot(filtered_data, "n_vertical", "Vertical Structure")


ggarrange(p_IVI, p_rich, p_stem_dens, p_vert, nrow = 2, ncol = 2, common.legend = T )


# filetr only final date
# Filter the data once
filtered_data_end <- df_indi_merged %>% 
  dplyr::filter(clim_scenario == "HISTO") %>% 
  dplyr::filter(year == 0 | year == 30)

# Create plots
p_IVI <- create_plot(filtered_data_end, "rIVI", "rIVI")
p_rich <- create_plot(filtered_data_end, "richness", "Richness")
p_stem_dens <- create_plot(filtered_data_end, "stem_density", "Stem Density")
p_vert <- create_plot(filtered_data_end, "n_vertical", "Vertical Structure")


ggarrange(p_IVI, p_rich, p_stem_dens, p_vert, nrow = 2, ncol = 2, common.legend = T )



# get Kruska test
# Perform Kruskal-Wallis Test for each indicator
kruskal_test_results <- filtered_data_end %>% 
 # filter(clim_scenario == "HISTO") %>%
  group_by(source) %>%
  summarise(
    rIVI_kruskal = kruskal.test(rIVI ~ clim_class)$p.value,
    richness_kruskal = kruskal.test(richness ~ clim_class)$p.value,
    stem_density_kruskal = kruskal.test(stem_density ~ clim_class)$p.value,
    n_vertical_kruskal = kruskal.test(n_vertical ~ clim_class)$p.value
  )



# Save results and tables -------------------------------------------------

# write field for RMarkdown
fwrite(df_field_ind_sub,  'outTable/df_field_indicators_landscapes.csv')
fwrite(df_sim_indicators, 'outTable/df_simulated_indicators.csv')
fwrite(df_sim_class,      'outTable/fin_sim_data.csv')
fwrite(df_compare_end,    'outTable/compare_field_sim_12_lands_end.csv')
fwrite(df_compare0,       'outTable/compare_field_sim_12_lands_start.csv')
fwrite(df_str_clustering_out,  'outTable/df_str_clustering_out_sankey.csv')

