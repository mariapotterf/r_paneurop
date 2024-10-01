
# Process: 
# Get indicators from simulated data
# - read iLand simulated data - simulated by Kilian on 29/06/2024
# - get vertical classes
# - calculate indicators form simulated data
# - average values between no seed and seed scenario - complete later by sensitivity analysis
# - evaluate they cchange change over time, how thtable they are given the clim cluster


# pre-analysis
# - inspect the development of 12 iLand landscapes


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
#library(cluster)

# Process input data -----------------------------------------------------------
# get simulated data
df_sim <- fread('outTable/df_simulated.csv')
head(df_sim)

# get field data to compare with simulated ones - stem density
df_field     <- fread('outData/veg_density_DBH.csv')
df_indicators <- fread('outData/indicators_for_cluster_analysis.csv') # summarized on plot level


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

# Create the data frame with the given pairs
df_sites_clusters <- data.frame(
  site = c("23_132", 
           "26_134", "15_133", "17_104", "22_101", "12_151",
           "24_146", "20_116", "12_117", "11_145", "19_160", "25_150"),
  cluster = c("1_1", "1_2", "1_3", "1_4", "2_1", "2_2",
              "2_3", "2_4", "2_5", "3_1", "3_2", "3_3")
)


# name the clusters from Kilian: if cluster analysis is run separately, the numbers do not fit! therefore,
# check the naming from plots (boxplot) an rename them to fit
df_sim_class <- df_sim %>%
  rename(landscape = cluster) %>% 
  mutate(clim_cluster = str_sub(landscape, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(landscape, -1, -1))  %>% # add indication of the strutural cluster (1,2,3,4,5) 
  mutate(clim_class = case_when(
    clim_cluster == 1 ~ "wet-warm-clay",  # Cluster 1: wet, cold, clay
    clim_cluster == 2 ~ "hot-dry-clay",  # Cluster 2: hot, dry, clay (more sand, less clay, more av.nitro than cluster 3)
    clim_cluster == 3 ~ "hot-dry-sand"    # Cluster 3: hot, dry, more sand
  ))  %>% 
  mutate(landscape_run = paste(landscape, run_nr))  # yunique run per lanscape scenario

str(df_sim_class)
head(df_sim_class)

df_sim_class <- df_sim_class %>% 
  mutate(unique_sim_run = paste(clim_model, clim_scenario, ext_seed, landscape_run, sep = "_"))  # yunique run per lanscape scenario




## filter field data for selected clusters (simulated landscapes) --------------
df_field_sub <- df_field %>%
  rename(landscape = cluster) %>%
  right_join(df_sites_clusters) %>%
  mutate(clim_cluster = str_sub(landscape, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(landscape, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)


# get Landsscape indicators for all of locations:
df_indicators <- df_indicators %>% 
 # rename(site = cluster) %>% 
  left_join(df_sites_clusters) %>%
 # mutate()
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1),
         clim_cluster_spei3 = factor(clim_cluster_spei3))   # add indication of the strutural cluster (1,2,3,4,5)

# make sure to interpret Kilinas clusters properly - if they fit with mine! 
# instect on field ata, using the tmp, prcp and clay content
df_field_ind_sub <- df_indicators %>% 
  dplyr::select(site,  rIVI, richness, stem_density, n_vertical, tmp, prcp, spei3, clay_extract, sand_extract,
                clim_cluster_spei3) %>% 
#  rename(site = cluster) %>% 
  right_join(df_sites_clusters) %>% 
  mutate(clim_cluster_Kilian = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster_Kilian = str_sub(cluster, -1, -1),
         clim_cluster_spei3 = factor(clim_cluster_spei3))  # add indication of the strutural cluster (1,2,3,4,5)



### Check clim clustering naming :-----------------------------------------------------------------------------

# Fiel data :
# Make boxplot to adjust naming of Kilians clusters:
p1 <- df_indicators %>% 
  ggplot(aes(x = clim_cluster_spei3,
             y = tmp)) +
  geom_boxplot() + 
  ggtitle("clim_clusters Maya (full)")

p2 <- df_indicators %>% 
  ggplot(aes(x = clim_cluster_spei3,
             y = prcp)) +
  geom_boxplot()

p3 <- df_indicators %>% 
  ggplot(aes(x = clim_cluster_spei3,
             y = spei3)) +
  geom_boxplot()

p4 <- df_indicators %>% 
  ggplot(aes(x = clim_cluster_spei3,
             y = clay_extract)) +
  geom_boxplot()

p5 <- df_indicators %>% 
  ggplot(aes(x = clim_cluster_spei3,
             y = sand_extract)) +
  geom_boxplot()

windows()
ggarrange(p1,p2,p3, p4, p5)





# Make boxplot to adjust naming of Kilians clusters:
p1 <- df_field_ind_sub %>% 
  ggplot(aes(x = clim_cluster_Kilian,
             y = tmp)) +
  geom_boxplot() +
  ggtitle("clim_clusters kilian (subset)")

p2 <- df_field_ind_sub %>% 
  ggplot(aes(x = clim_cluster_Kilian,
             y = prcp)) +
  geom_boxplot()

p3 <- df_field_ind_sub %>% 
  ggplot(aes(x = clim_cluster_Kilian,
             y = spei3)) +
  geom_boxplot()

p4 <- df_field_ind_sub %>% 
  ggplot(aes(x = clim_cluster_Kilian,
             y = clay_extract)) +
  geom_boxplot()

p5 <- df_field_ind_sub %>% 
  ggplot(aes(x = clim_cluster_Kilian,
             y = sand_extract)) +
  geom_boxplot()

windows()
ggarrange(p1,p2,p3, p4, p5)




# Process simulated data --------------------------------------------------------
## filter teh data for teh baseic scenario: one lansca, one clim scenarions, ... --------------------
# the most basic example:

df_sim %>% 
  dplyr::filter(clim_model == 'ICHEC' & clim_scenario == 'HISTO' &
                  cluster == '1_4' & 
                  run_nr == '1'&
                  ext_seed  == 'noseed' &
                  count_ha > 0
                ) %>% 
  #View()#
  ggplot(aes(x = year,
             y = count_ha,
             color = interaction(species, category  ),
             group = interaction(species, category  ))) +
  geom_line() +
    facet_wrap(species~category)




## Inspect clim drivers: For all sites  -------------------
# to understand which cluster is what:
# Create an empty list to store the plots
plot_list_all <- list()

# Loop through combinations and create plots
for (combination in plot_combinations) {
  p <- ggplot(df_indicators, aes_string(x = combination$x, y = combination$y)) +
    geom_point(color = 'grey', alpha = 0.3) + 
   # geom_point(data = df_field_ind_sub, aes_string(x = combination$x, y = combination$y, color = clim_cluster), size = 4) +
   # geom_smooth(method = "loess", se = FALSE) +
   # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    labs(title = paste(combination$x, "vs", combination$y),
         x = combination$x,
         y = combination$y) +
    theme_bw() +
    theme(aspect.ratio = 1)  # Ensure the plot is square
  
  plot_list_all[[paste(combination$x, combination$y, sep = "_vs_")]] <- p
}
# Arrange and print all plots using ggarrange
ggarrange(plotlist = plot_list_all, ncol = 2, nrow = 3, common.legend = TRUE) # 
 

# make a plot manually: tmp spei3
p1 <- ggplot(df_indicators, aes(x = tmp, y = spei33)) +
  geom_point(color = 'grey', alpha = 0.3) + 
  geom_text(data = df_field_ind_sub, aes(x = tmp, y = spei3, label = cluster, color = clim_cluster), 
            vjust = -0.1, size = 6) +
  geom_point(data = df_field_ind_sub, aes(x = tmp, y = spei3, color = clim_cluster), size = 3) +
  # labs(#title = paste(combination$x, "vs", combination$y),
  #   x = "tmp",
  #   y = "spei3") +
  theme_bw() +
  theme(aspect.ratio = 1)  # Ensure the plot is square


p2 <- ggplot(df_indicators, aes(x = tmp_z, y = spei3)) +
  geom_point(color = 'grey', alpha = 0.3) + 
  geom_text(data = df_field_ind_sub, aes(x = tmp_z, y = spei3, label = cluster, color = clim_cluster), 
            vjust = -0.1, size = 6) +
  geom_point(data = df_field_ind_sub, aes(x = tmp_z, y = spei3, color = clim_cluster), size = 3) +
   theme_bw() +
  theme(aspect.ratio = 1)  # Ensure the plot is square


p3 <- ggplot(df_indicators, aes(x = tmp_z, y = prcp_z)) +
  geom_point(color = 'grey', alpha = 0.3) + 
  geom_text(data = df_field_ind_sub, aes(x = tmp_z, y = prcp_z, label = cluster, color = clim_cluster), 
            vjust = -0.1, size = 6) +
  geom_point(data = df_field_ind_sub, aes(x = tmp_z, y = prcp_z, color = clim_cluster), size = 3) +
  theme_bw() +
  theme(aspect.ratio = 1)  # Ensure the plot is square

ggarrange(p1,p2,p3, common.legend = T, ncol = 3, nrow = 1)

# plot only subset: climatic clusters (12 points)
df_field_ind_sub %>% 
  ggplot(aes(x = spei3,
                 y = tmp_z)) +
  geom_point(aes(color = clim_cluster)) + 
  geom_smooth()


# the clim cluster 2_1  seems a bit strange: get the characteristic of clusters 

# # inspect single example
# df_sub1 <- df_sim %>% 
#   dplyr::filter(cluster == "1_1") %>% 
#   dplyr::filter(clim_model == "NCC" & clim_scenario == "HISTO" & ext_seed == "noseed")
# 
# head(df_sub1)

## Inspect simulated data: ------------------------------------------------------
# 12 clusters, 
# 3 climatic models
# 4 climate scenarios
# 2 seed added/no
unique(df_sim2$landscape)  # 12, 3 climatic, 4 strcutural in each climate cluster
#[1] "1_1" "1_2" "1_3" "1_4" "2_1" "2_2" "2_3" "2_4" "2_5" "3_1" "3_2" "3_3"

unique(df_sim2$species) # 24
#[1] "acps" "cabe" "potr" "saca" "quro" "acca" "frex" "ulgl" "lade" "tico" "piab" "fasy" "psme"
#[14] "soau" "pisy" "bepe" "abal" "alin" "algl" "soar" "coav" "acpl" "rops" "casa"


unique(df_sim2$clim_model)
# [1] "ICHEC" "MPI"   "NCC"  

unique(df_sim2$clim_scenario)
# [1] "HISTO" "RCP26" "RCP45" "RCP85"

# external seed: external seed present or not?? TRUE/FALSE
# years: 1-30
unique(df_sim2$run_nr)  # 5 repetitions



### Indicators from simulated ata : ------------------------------------------------------ 

# richness
# rIVI
# stem density
# vertical classes
# summarize across simulation repetition and across the clim model

##### Investigate on average seed scenario -----------------------------------------------

# Calculate averages grouped by year, seed strategy, clim_cluster, species, run_nr, and vertical category
df_avg <- df_sim_class %>%
  group_by(year, landscape_run, clim_cluster, clim_class, species, category) %>%
  summarise(mean_count_ha = mean(count_ha, na.rm = TRUE),
            mean_basal_area = mean(basal_area_m2, na.rm = TRUE))


# START - average seed scenarios 

###### Species richness ------------------------------------
df_richness <- 
  df_avg %>% 
  group_by(year, landscape_run, clim_cluster, clim_class) %>%  #clim_modelclim_cluster, 
  dplyr::filter(mean_count_ha>0) %>% 
  summarize(richness = n_distinct(species))# %>% 
  #View()


df_richness %>% 
   dplyr::filter(landscape_run == "1_1 5" )  %>%
  #View()
  ggplot() + 
  geom_line(aes(x = year,
                y = richness,
                color = landscape_run,
                group = landscape_run)) + 
  facet_grid(.~ clim_cluster)


###### Species importance value (relative, from rel_density and rel_BA) ------------------------------------
# get species importance value: from relative density, relative BA
# first calculate the total values per ha/landscape/run, then add it to original table to calculate teh rIVI based on relative dominance
df_sum_landscape <- df_avg %>% 
  group_by(year,  landscape_run, clim_cluster, clim_class) %>%# clim_cluster clim_model ,, clim_scenario, ext_seed,
  summarize(sum_stems = sum(mean_count_ha, na.rm = T),
            sum_BA = sum(mean_basal_area , na.rm = T)) %>% 
  ungroup()


df_IVI <- df_avg %>% 
  # group by speies across the levels
  group_by(year, species, landscape_run, clim_cluster, clim_class) %>% #clim_cluster,clim_model , landscape, clim_scenario, ext_seed
  summarize(sp_dens = sum(mean_count_ha),
            sp_BA   = sum(mean_basal_area)) %>% 
  ungroup(.) %>% 
  left_join(df_sum_landscape) %>% # merge sums per whole cluster
  mutate(rel_dens  = sp_dens/sum_stems,
         rel_BA  = sp_BA/sum_BA,
         rIVI = ( rel_dens +rel_BA)/2) %>% # relative species importance value
  # filter teh dominant species - highest rIVI
  group_by(year, landscape_run) %>% #clim_model ,,  clim_scenario, ext_seed
  dplyr::filter(rIVI == max(rIVI)) %>%
  sample_n(1) %>%  # Select a random row if there are ties
  rename(dominant_species      = species)




##### Structure: Vertical classes & stem density -------------------------------------------------------------
# inspeact vertical classes
df_structure <- df_avg %>% 
  group_by(year, landscape_run, clim_cluster, clim_class ) %>% #clim_model , , clim_scenario, ext_seed
  summarize(n_vertical = n_distinct(category),
            stem_density = sum(mean_count_ha)) 


# change in number of vertical layers
ggplot(df_structure) + 
  geom_line(aes(x = year,
                y = n_vertical,
                color = clim_cluster ,
                group = landscape_run
                ),
            alpha = 0.5) #+ 
#  facet_wrap(clim_scenario~ext_seed)



# change in stem density
ggplot(df_structure) + 
  geom_line(aes(x = year,
                y = stem_density,
                color = clim_cluster ,
                group = landscape_run
            ), alpha = 0.5 ) +
  
  facet_wrap(clim_cluster~.)


# merge df indicators 
df_sim_indicators <-
  df_richness %>% 
  left_join(df_IVI, by = join_by(year, landscape_run, clim_cluster, clim_class)) %>% #, clim_scenario, ext_seed
  left_join(df_structure, by = join_by(year, landscape_run, clim_cluster, clim_class)) %>% #, clim_scenario, ext_seed
#  mutate(clim_cluster = str_sub(landscape, 1, 1),  # add indication of the climatic landscape (1,2,3)
#         str_cluster = str_sub(landscape, -1, -1))  #%>% # add indication of the strutural cluster (1,2,3,4)

#  %>%   #%>% # add indication of the strutural cluster (1,2,3,4)
  mutate(site = str_sub(landscape_run, 1,3))  # remove 'run' indication


head(df_sim_indicators)








##### All simulations run  -------------------------
###### Species richness ------------------------------------

# Get grouping table : 1380 unique simulations combination
df_simulations_groups <- df_sim_class %>% 
  dplyr::select(year, clim_model,clim_cluster,str_cluster, ext_seed,landscape_run, unique_sim_run) %>% 
  distinct() %>%  # remove duplicated rows
  mutate(landscape = paste(clim_cluster,str_cluster, sep = "_"))
# nrows with years: 42720


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
  left_join(df_structure) #%>% 
  #mutate(clim_cluster_kilian = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
   #      str_cluster_Kilian = str_sub(cluster, -1, -1)) 

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




##### evaluate initial state with my sites (only 12 sites! )   -------------------------
# filterr initial state: year == 0
df_sim_indicators0 <- df_sim_indicators %>% 
  ungroup() %>% 
  dplyr::filter(year == 0 ) %>% #& clim_scenario == "HISTO" %>% 
  dplyr::select(landscape_run,  rIVI, richness, stem_density, n_vertical, landscape) #%>% 
  


###### merge field data with simulated data in year 0 
df_compare <- 
  df_field_ind_sub %>% 
  left_join(df_sim_indicators0, by = c( "cluster" = "landscape"), suffix = c("_field", "_simul")) %>% 
  dplyr::select(cluster, site, ends_with("_field"), ends_with("_simul")) %>% 
 # mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
 #        str_cluster = str_sub(cluster, -1, -1))  %>% # add indication of the strutural cluster (1,2,3,4)
 na.omit() # remove empty one
  
  

print(df_compare)


# understand structural clusters :  ------------------
p1 <- df_field_ind_sub %>% 
  ggplot(aes(x = cluster,
             y = rIVI,
             fill = clim_cluster_spei3)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- df_field_ind_sub %>% 
  ggplot(aes(x = cluster,
             y = richness,
             fill = clim_cluster_spei3)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <-df_field_ind_sub %>% 
  ggplot(aes(x = cluster,
             y = stem_density,
             fill = clim_cluster_spei3)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p4<-df_field_ind_sub %>% 
  ggplot(aes(x = cluster,
             y = n_vertical,
             fill = clim_cluster_spei3)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggarrange(p1,p2,p3,p4, common.legend = T)


# explain climate and soil
p.tmp <- df_field_ind_sub %>% 
  ggplot(aes(x = clim_cluster_spei3,
             y = tmp)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))

p.prcp <- df_field_ind_sub %>% 
  ggplot(aes(x = clim_cluster_spei3,
             y = prcp)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))

p.spei <- df_field_ind_sub %>% 
  ggplot(aes(x = clim_cluster_spei3,
             y = spei3)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))


ggarrange(p.tmp, p.prcp, p.spei,common.legend = T)

# check indicators
df_field_ind_sub %>% 
  dplyr::select(stem_density, n_vertical, richness, rIVI, cluster, prcp, tmp ) %>% 
  arrange(stem_density)



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
    composite_score = round(scaled_stem_density + scaled_n_vertical + scaled_richness - scaled_rIVI, 2),
    cluster = factor(cluster, levels = cluster[order(composite_score, decreasing = TRUE)])
  ) %>% 
  mutate(composite_class = factor(composite_score)) %>% 
  mutate(clim_class = case_when(
    clim_cluster_spei3  == 1 ~ "wet-warm-clay",  # Cluster 1: wet, cold, clay
    clim_cluster_spei3  == 2 ~ "hot-dry-clay",  # Cluster 2: hot, dry, clay (more sand, less clay, more av.nitro than cluster 3)
    clim_cluster_spei3  == 3 ~ "hot-dry-sand"    # Cluster 3: hot, dry, more sand
  ))  #%>% 
  

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
  dplyr::select(site,  rIVI, richness, stem_density, n_vertical, clim_cluster_spei3, clim_class) %>%
  mutate(year = 1) %>% 
  mutate(type = 'field',
         landscape_run = 'field') %>% 
  rename(clim_cluster = clim_cluster_spei3)


df_simulated <- df_sim_indicators %>% 
  ungroup() %>%  
  mutate(type = 'simulated') %>% 
  dplyr::select(landscape,  rIVI, richness, stem_density, n_vertical, clim_cluster,   clim_class,  year,   type, landscape_run)
 

head(df_indicators_field)

head(df_simulated)


df_merge <- rbind(df_simulated, df_indicators_field)


df_merge %>% 
  ggplot(aes(x = year,
             y = stem_density,
             group = landscape_run)) + 
  geom_line()



# TEST START - for one variable

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
#TEST END




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
fwrite(df_field_ind_sub, 'outTable/df_field_indicators_landscapes.csv')
fwrite(df_sim_indicators, 'outTable/df_simulated_indicators.csv')
fwrite(df_sim_class, 'outTable/fin_sim_data.csv')

