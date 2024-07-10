

# read iLand simulated data
# simulated by Kilian on 29/06/2024
# calculate my variables
# evaluate how my indicators will change over time
# classify teh vertical classes


# check data:
# inspect my variables with the simulated ones in year 0/1

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

df <- fread('outTable/df_simulated.csv')
View(df)


# inspect single example
df_sub1 <- df_cl %>% 
  dplyr::filter(cluster == "1_1") %>% 
  dplyr::filter(clim_model == "NCC" & clim_scenario == "HISTO" & ext_seed == "FALSE")

View(df_sub1)


# Check simulation data:
# 12 clusters, 
# 3 climatic models
# 4 climate scenarios
# 2 seed added/no
unique(df$cluster)  # 12, 3 climatic, 4 strcutural in each climate cluster
#[1] "1_1" "1_2" "1_3" "1_4" "2_1" "2_2" "2_3" "2_4" "2_5" "3_1" "3_2" "3_3"

unique(df$species) # 24
#[1] "acps" "cabe" "potr" "saca" "quro" "acca" "frex" "ulgl" "lade" "tico" "piab" "fasy" "psme"
#[14] "soau" "pisy" "bepe" "abal" "alin" "algl" "soar" "coav" "acpl" "rops" "casa"


unique(df$clim_model)
# [1] "ICHEC" "MPI"   "NCC"  

unique(df$clim_scenario)
# [1] "HISTO" "RCP26" "RCP45" "RCP85"

# external seed: external seed present or not?? TRUE/FALSE
# years: 1-30
unique(df$run_nr)  # 5 repetitions


# need to run vertical classification:
# mature > 10 cm dbh
# saplings 20cm -2 m
# juveniles > 2m - 10 cm dbh
# calculate my summary variables
df <- df %>%
  mutate(clim_cluster = str_sub(cluster, 1, 1))  # add indication of the climatic cluster (1,2,3)

# df_cl %>% 
#   dplyr::filter(category == "unclassified") %>% 
#   summary()

# list of teh sites vs Kilian's clusters: 
# #Site - Cluster
# 
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


# calculate my variables on site level: ------------------------------------------------------ 
# richness
# rIVI
# stem density
# vertical classes

### Species richness ------------------------------------
df_richness <- df %>% 
  group_by(year, clim_cluster, cluster, clim_model , clim_scenario, ext_seed, run_nr) %>% 
  summarize(richness = n_distinct(species))

df_richness %>% 
  dplyr::filter(run_nr == 5)  %>% 
  ggplot() + 
  geom_line(aes(x = year,
                y = richness,
                color = cluster,
                group = cluster)) + 
  facet_wrap(clim_scenario~ext_seed)


### Species importance value (relative, from rel_density and rel_BA) ------------------------------------
# get species importance value: from relative density, relative BA
# first calculate the total values per ha, then add it to original table to calculate teh rIVI based on relative dominance
df_sum_cluster <- df_cl %>% 
  group_by(year, clim_cluster, cluster, clim_model , clim_scenario, external_seed) %>% 
  summarize(sum_stems = sum(count_ha, na.rm = T),
            sum_BA = sum(basal_area_m2, na.rm = T)) %>% 
  ungroup()


df_IVI <- df_cl %>% 
  # group by speies across the levels
  group_by(year, species,clim_cluster, cluster, clim_model , clim_scenario, external_seed) %>% 
  summarize(sp_dens = sum(count_ha),
            sp_BA = sum(basal_area_m2)) %>% 
  ungroup(.) %>% 
    left_join(df_sum_cluster) %>% # merge sums per whole cluster
  mutate(rel_dens  = sp_dens/sum_stems,
         rel_BA  = sp_BA/sum_BA,
         rIVI = ( rel_dens +rel_BA)/2) %>% # relative species importance value
  # filter teh dominant species - highest rIVI
  group_by(year, cluster, clim_model , clim_scenario, external_seed) %>% 
  dplyr::filter(rIVI == max(rIVI)) %>%
  sample_n(1) #%>%  # Select a random row if there are ties
  
  
# inspect values by filtering: seems ok
df_IVI %>% 
  dplyr::filter(cluster == '1_4' & clim_scenario == "HISTO" & clim_model == "NCC" & external_seed == "FALSE") %>% 
  View()
  
  
  

# Structure: Vertical classes & stem density -------------------------------------------------------------
# inspeact vertical classes
df_structure <- df_cl %>% 
  group_by(year, cluster, clim_model , clim_scenario, external_seed) %>% 
  summarize(n_vertical = n_distinct(category),
            stem_dens = sum(count_ha)) 

# change in number of vertical layers
ggplot(df_structure) + 
  geom_line(aes(x = year,
                y = n_vertical,
                color = cluster,
                group = cluster)) + 
  facet_wrap(clim_scenario~external_seed)



# change in stem density
ggplot(df_structure) + 
  geom_line(aes(x = year,
                y = stem_dens,
                color = cluster,
                group = cluster)) + 
  facet_wrap(clim_scenario~external_seed)


  