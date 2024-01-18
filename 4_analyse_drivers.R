# create final database
# merge veg data
# predictors - climate, disturbance chars, soils
# dependent - regeneretion, density, richness
# output file
# 

# Libs --------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(terra)
library(ggplot2)
library(lubridate)
library(stringr)


library(dunn.test) # Duncan test, post hoc Kruskal-Wallis


# stats 
library(MASS) #glm.nb for negative binomial models
library(glmmTMB) #Salamanders dataset and lots of families


library(lme4) #pseudo-R2 and and mixed model functionality
library(MuMIn) #dredge function for comparing models, AIC, marginal R^2 for mixed models
library(sjmisc) #pseudo R^2 - but I think gets loaded as a dependency
library(DHARMa) #model diagnostics
library(effects) #what do my marginal effects look like?
library(performance) #binomial model diagnostics
library(emmeans) #post hoc for categorical predictors


# read data --------------------------------------------------------------------

climate           <- fread("outData/climate_18_23.csv")
spei              <- fread("outData/SPEI.csv")
distance_to_edge  <- fread("outData/distance_to_edge.csv")
disturbance_chars <- fread("outData/disturbance_chars.csv")
soil              <- as.data.frame(terra::vect("rawData/extracted_soil_data/extracted_soil_data.gpkg"))
terrain           <- fread("outData/terrain.csv")

# get vegetation data
load("outData/veg.Rdata")

# clean up data -----------------------------------------------------------------
# sum stems first up across vertical groups!!
df_stems <- stem_dens_ha_cluster %>% 
  group_by(cluster, new_manag) %>% 
  summarise(sum_stems = sum(total_stems_all_species))# %>% 


disturbance_chars <- disturbance_chars %>%
  mutate(country = as.character(country))

# merge all predictors together, ID level
df_predictors <- 
  climate %>% 
  left_join(select(distance_to_edge, c(-country, -dist_ID ))) %>% 
  left_join(select(disturbance_chars, c(-region, -country ))) %>% 
  left_join(soil) %>%
  left_join(select(terrain, c(-country))) 


# select the dominant species per cluster
df_IVI_max <- df_IVI %>% 
  dplyr::select(cluster, Species, new_manag, country,rIVI) %>% 
  ungroup(.) %>% 
  group_by(cluster, country, new_manag) %>% 
  filter(rIVI == max(rIVI)) %>% 
  slice(1)

# cluster level
df_cluster <-
  df_IVI %>% 
  left_join(df_richness) %>%
  left_join(df_vert) %>%
  left_join(df_stems)# %>%


# Analysis ------------------------------------------------------------
# 1. How many clusters are only manag/unmanag/mixed?
manag_types <- stem_dens_species_long %>%
  ungroup(.) %>% 
  dplyr::select(cluster, new_manag) %>% 
  group_by(cluster, new_manag) %>% 
  distinct(.) %>%
  ungroup(.) %>% 
  group_by(cluster) %>% 
  summarise(type = case_when(
    all(new_manag == 'Unmanaged') ~ 'Unmanaged',
    all(new_manag == 'Managed') ~ 'Managed',
    TRUE ~ 'Mixed'
  )) %>%
  ungroup(.) %>% 
  dplyr::select(-cluster) %>% 
  table()




#MAN MIX UNM 
#616 149  84 

## 2. species composition per stems: saplings vs juveniles?? --------------
# how many plots with Mature trees?
# 
dom_species <- stem_dens_species_long %>% 
  filter(VegType != 'Survivor' ) %>%  # Survivors: on 191
  group_by(Species, VegType, new_manag ) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup(.) %>% 
  group_by(VegType) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share = sum_stems/sum_vegType*100) #%>% 
  

# see all species
dom_species %>% 
 
  #dplyr::filter(share > 1) %>% 
  ggplot(aes(x = Species,
             y = share,
             fill = VegType)) +
  geom_bar(position = "fill", stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
  
  #geom_col('')



# all species by aboundance
dom_species %>% 
  #dplyr::filter(share > 1) %>% 
  ggplot(aes(x = Species,
             y = share,
             fill = VegType)) +
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#geom_col('')



# species by management
dom_species %>% 
  #dplyr::filter(VegType == 'Regeneration') %>% 
  dplyr::filter(share > 0) %>% 
  ggplot(aes(x = reorder(Species, -share),
             y = share,
             fill = VegType)) +
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(VegType~new_manag)
#geom_col('')


# 3. how many clusters/plots has less then 2000 regen/ha? ----------------------------------
# which country is the most affected? (where they ccur most often?)



df_stock_density <- 
  df_stems %>% 
  mutate(dens_category = case_when(sum_stems < 500 ~ '<500',
                                   sum_stems >= 500 &sum_stems < 1000 ~ '500-1000',
                                   sum_stems >= 1000 &sum_stems < 1500 ~ '1000-1500',
                                   sum_stems >= 1500 &sum_stems < 2000 ~ '1500-2000',
                                   sum_stems >= 2000 &sum_stems < 2500 ~ '2000-2500',
                                   sum_stems >= 2500 ~ '>2500'
                                   )) %>% 
  mutate(dens_category = factor(dens_category, levels = 
                                  c('<500','500-1000','1000-1500', '1500-2000',  '2000-2500', '>2500' )))
str(df_stock_density)

df_stock_density %>% 
  ggplot(aes(x = new_manag,
             fill = dens_category)) +
  geom_bar(position = 'fill') +
  ylab('Share [%]') + 
  theme_classic()
  

# Compare categories occurences by chi square


### 4. stems vs weather - mean 2019-2022, or SPEI ---------------------------------
# average predictors by cluster
df_predictors_clust <- df_predictors %>% 
  mutate(cluster = str_sub(ID, 4, -3)) %>% 
  group_by(cluster) %>%
  filter(year %in% 2021:2023) %>% 
  summarise(tmp                  = mean(tmp, na.rm = T),
            prec                 = mean(prec, na.rm = T), 
            tmp_z                = mean(tmp_z, na.rm = T),
            prcp_z               = mean(prcp_z, na.rm = T),
            distance             = mean(distance, na.rm = T),
            disturbance_year     = mean(disturbance_year, na.rm = T),
            disturbance_severity = mean(disturbance_severity, na.rm = T),
            disturbance_agent    = mean(disturbance_agent, na.rm = T),
            elevation            = mean(elevation, na.rm = T),
            sand_extract         = mean(sand_extract, na.rm = T),
            silt_extract         = mean(silt_extract, na.rm = T),
            clay_extract         = mean(clay_extract, na.rm = T),
            depth_extract= mean(depth_extract, na.rm = T),
            av.nitro  = mean(av.nitro, na.rm = T),
            slope = mean(slope, na.rm = T),
            aspect= mean(aspect, na.rm = T))
  
length(unique(df_predictors_clust$cluster))

  

### get average spei after disturbances ------------------------------------------
# need to update, as too many points are missing!!!
# eg ID 17_21_114_1

#spei_sub <- 
  spei %>% 
  dplyr::mutate(year  = lubridate::year(date),
                month = lubridate::month(date)) %>%  
  dplyr::select(-c(date, -scale)) %>% 
  mutate(cluster = str_sub(ID, 4, -3)) %>%
    ungroup(.) %>% 
    dplyr::filter(year %in% 2021:2023) %>%
  View()
    group_by(cluster) %>%
    summarise(spei = mean(spei, na.rm = T)) #%>% 
    
    
# join to stems counts: filter, as from veg data I have already excluded Italy and teh uncomplete clusters
df_fin <- df_stems %>%
  left_join(df_predictors_clust, by = join_by(cluster)) %>% 
  #left_join(spei_sub, by = join_by(cluster)) %>% 
  mutate(new_manag = as.factor(new_manag)) %>% 
  mutate(sum_stems = round(sum_stems)) %>% 


length(unique(df_fin$cluster))


pairs(sum_stems   ~ tmp + tmp_z + distance, df_fin)

# some values are missing predictors: eg. 18_178, 18_184

#5. test model ---------------------------------

# test glm with neg bin family
library(MASS)
m.negbin <- glm.nb(sum_stems ~ tmp  + tmp_z +  distance + new_manag, data = df_fin, na.action = "na.fail") # + tmp_z 

dredge(m.negbin)

# zero inflated model
library(pscl)
m.zip <- zeroinfl(sum_stems ~ tmp + tmp_z + distance + new_manag, 
                  dist = "poisson", data = df_fin)


# spei has infinity values?? eg here: 21_142


# 6. management effect: compare manag vs unmanaged stem densities, richness

ggplot(data = df_richness,
       aes(x = new_manag, y = richness)) + 
  # Add bars as medians
  stat_summary(fun = "median", 
               geom = "bar", 
               alpha = .7) +
  stat_summary(
    data = df_richness,
    mapping = aes(x = new_manag, y = richness),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) #+
  #coord_cartesian(ylim=c(60, 67)) # ylim=c(59,66)


# compre densities 
stem_dens_ha_cluster %>% 
  group_by(new_manag) %>% 
  summarise(mean = mean(total_stems_all_species),
            sd = sd(total_stems_all_species),
            median = median(total_stems_all_species))


windows()
ggplot(df_stems, 
       aes(x = new_manag, y = sum_stems )) + 
  # Add bars as medians
  stat_summary(fun = "median", 
               geom = "bar", 
               alpha = .7) +
  stat_summary(
    data = df_stems,
    mapping = aes(x = new_manag, y = sum_stems),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) +
  ylab('Stems [ha]') + 
  theme_classic()
#coord_cartesian(ylim=c(60, 67)) # ylim=c(59,66)

# test significance - for 2 groups only
wilcox.test(sum_stems ~ new_manag, data = df_stems)


# for several groups:
kruskal.test(sum_stems ~ new_manag, data = df_stems)


# Perform the Dunn test
dunnTest <- dunn.test(df_stems$sum_stems, df_stems$new_manag, method="bonferroni")
dunnTest
  
            