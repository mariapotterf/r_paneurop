

# extract variables 
# density
# tree species
# IVI
# how to summarize data? average oevr clusters???
# smart way to keep empty clusters?
# make a major df


# analyse structure and composition
# recalculate density to hectares
# correctly account for the empty plots!! get a master df with all 


# what predictors to get? 
# IVI - species importance value, based on stem density and basal area
# # include management


# what is 'cluster?  = region +group?
# country = '
# management = fill in missing values correctly: based on the first estimation per plot (without completed species)
gc()

library(data.table)
library(dplyr)
library(ggplot2)
library(terra)
library(tidyr)


dat_raw <- fread('rawData/working_directory/rapid_assessment_mdf.csv')

# check: has 4 mature trees? in ID 12_17_112_3

#dat %>% 
 # dplyr::filter(ID == "12_17_112_3") %>% View()

# remove Italy and Vaia - windthrow
dat <- dat_raw %>% 
  dplyr::filter( country != 17) #   ~ "IT",  # Italy)

# create management based on Christians scrips;
# get management intensity; for each plot and for each cluster: per plot, max intensity is 5, per cluster it is 25 (5*5)
dat2 <- dat %>% 
  dplyr::filter(dist == TRUE) %>% # remove the plot if not disturbed
  mutate(cluster = paste(region, group, sep = '_')) %>% 
  # recode teh variables as tehy have numbers 1:3, not truie/false
  mutate(planting = case_when(
    planting == 1 ~ NA,       # Recode 1 to NA
    planting == 2 ~ TRUE,     # Recode 2 to TRUE
    planting == 3 ~ FALSE,    # Recode 3 to FALSE
    is.na(planting) ~ NA,     # Preserve existing NA values - converrt to 0 later
    TRUE ~ planting            # Keep other values as they are
  )) %>% 
  mutate(anti_browsing = case_when(
    anti_browsing == 1 ~ NA,       # Recode 1 to NA
    anti_browsing == 2 ~ TRUE,     # Recode 2 to TRUE
    anti_browsing == 3 ~ FALSE,    # Recode 3 to FALSE
    is.na(anti_browsing) ~ NA,     # Preserve existing NA values
    TRUE ~ anti_browsing            # Keep other values as they are
   )) %>%
  mutate(windthrow = case_when(
    windthrow == 1 ~ NA,       # Recode 1 to NA
    windthrow == 2 ~ TRUE,     # Recode 2 to TRUE
    windthrow == 3 ~ FALSE,    # Recode 3 to FALSE
    is.na(windthrow) ~ NA,     # Preserve existing NA values - converrt to 0 later
    TRUE ~ windthrow            # Keep other values as they are
  )) %>% 
  mutate(deadwood = case_when(
    deadwood == 1 ~ NA,       # Recode 1 to NA
    deadwood == 2 ~ TRUE,     # Recode 2 to TRUE
    deadwood == 3 ~ FALSE,    # Recode 3 to FALSE
    is.na(deadwood) ~ NA,     # Preserve existing NA values - converrt to 0 later
    TRUE ~ deadwood            # Keep other values as they are
  )) %>% 
   mutate(
    logging_trail = ifelse(!is.na(logging_trail) & logging_trail == TRUE, 1, 0),
    clear         = ifelse(!is.na(clear) & clear == TRUE, 1, 0),
    grndwrk       = ifelse(!is.na(grndwrk) & grndwrk == TRUE, 1, 0),
    planting      = ifelse(!is.na(planting) & planting == TRUE, 1, 0),
    anti_browsing = ifelse(!is.na(anti_browsing) & anti_browsing == TRUE, 1, 0),
    manag_intensity      = logging_trail + clear + grndwrk + planting + anti_browsing,  # get management intensity: rate on cluster cluster level
    salvage_intensity    = logging_trail + clear + grndwrk,  # get salvage intensity: how much the site was altered by harvest?rate on cluster cluster level
    protection_intensity = planting + anti_browsing  # get were trees plantedor even fenced? rate on cluster cluster level
)


#View(dat_test)
dat2 <- dat2 %>% 
  mutate(country = case_when(
    country == 11 ~ "DE",  # Germany
    country == 12 ~ "PL",  # Poland
    country == 13 ~ "CZ",  # Czech Republic
    country == 14 ~ "AT",  # Austria
    country == 15 ~ "SK",  # Slovakia
    country == 16 ~ "SI",  # Slovenia
    country == 17 ~ "IT",  # Italy
    country == 18 ~ "CH",  # Switzerland
    country == 19 ~ "FR",  # France
    TRUE ~ NA_character_      # Fallback in case of unidentified country_id
  ))
# Name countries
#11 germany
#12 poland
#13 czech
#14 austria
#15 slovakia
#16 slovenia
#17 italy
#18 switzerland
#19 france


#change names of veg classes
dat2 <- dat2 %>% 
  mutate(VegType = case_when(
    VegType == "Regeneration" ~ "Saplings",  
    VegType == "advRegeneration" ~ "Juveniles",
    VegType == "Survivor" ~ "Mature",
    TRUE ~ NA_character_      # Fallback in case of unidentified country_id
  ))  

dat2 <- dat2 %>% 
  mutate(VegType  = factor(VegType, levels = c("Mature",
                                               "Juveniles",
                                               "Saplings"
                                               )))
  
# filetr:
# 13_26_165_3
# 13_26_165_5
# 13_26_165_1
# 13_26_165_2
# 13_26_165_4

d2_select_opocno <- dat2 %>% 
  dplyr::filter(ID %in% c("13_26_165_3",
                           "13_26_165_5",
                          "13_26_165_1",
                           "13_26_165_2",
                          "13_26_165_4")) %>% 
  dplyr::select(Species, VegType, Variable, value) %>% 
  dplyr::filter(Variable == "n") %>% 
  na.omit()

#View(d2_select_opocno)

# filter sub-plots that are misplaced/double recorded
# visually checked on May 8th, 2024
remove_sub_plots <- c('18_22_110_7' ,
                      '18_22_118_6', 
                      '13_15_145_4', 
                      '13_26_142_2', 
                      '11_14_123_2', 
                      '11_18_106_6', 
                      '11_18_178_5', 
                      '11_25_128_3', 
                      '12_17_107_4', 
                      '12_17_143_5' )

#   cluster       ID            desc
# "22_110" - ID: 18_22_110_7 - out of +
#   "22_118" - ID: 18_22_118_6 - out of +
#   "15_145" - ID: 13_15_145_4 - missing point geometry
# "26_142" - ID: 13_26_142_2 - overlying on 1
# "14_123" - ID: 11_14_123_2 - overly on 2
# "18_106" - ID: 11_18_106_6 - out of +
#   "18_178" - ID: 11_18_178_5 - out of +
#   "25_128" - ID: 11_25_128_3 - overlaying 2
# "17_107" - ID: 12_17_107_4 - out of +
#   "17_143" - ID: 12_17_143_5 - out of +

dat2 <- dat2 %>% 
  dplyr::filter(!ID %in% remove_sub_plots )
  
# check n_plots per cluster
#table(dat$)


# Create master df with empty plots - eg no trees found on them
df_master <- 
  dat2 %>% 
  dplyr::select(country,region, group, cluster, point, ID) %>%  #  region, group, 
  unique() %>%  # remove duplicated rows
  group_by(country, region, group, cluster) %>%  # region, group,
  summarize(n_plots = n())

length(unique(df_master$cluster))


# get cluster ID for the clusters with 6 plots = are they marked as disturbance yes?

cluster_6_plots <- df_master %>% 
  dplyr::filter(n_plots == 6) %>% 
  distinct(cluster) %>% 
  pull()

cluster_6_plots


# dat %>% 
#   dplyr::filter(cluster %in% cluster_6_plots) %>% 
#   View()
# 
# table(df_master$n_plots)
# 

# keep clusters with 4, as some disturbace plots were very small
# remove if there is less records
dat2 <- dat2 %>% 
  right_join(df_master, by = join_by(country, region, group, cluster)) %>%
  dplyr::filter(n_plots > 3)


length(unique(dat2$cluster)) # 849

## get individual management types ---------------------------
df_individual_management <- dat2 %>% 
  dplyr::filter(dist == TRUE) %>% # remove the plot if not disturbed
  mutate(cluster = paste(region, group, sep = '_')) %>% 
  mutate(
    logging_trail = ifelse(!is.na(logging_trail) & logging_trail == TRUE, 1, 0),
    clear         = ifelse(!is.na(clear) & clear == TRUE, 1, 0),
    grndwrk       = ifelse(!is.na(grndwrk) & grndwrk == TRUE, 1, 0),
    planting      = ifelse(!is.na(planting) & planting == TRUE, 1, 0),
    anti_browsing = ifelse(!is.na(anti_browsing) & anti_browsing == TRUE, 1, 0),
    windthrow = ifelse(!is.na(windthrow) & windthrow == TRUE, 1, 0),
    deadwood = ifelse(!is.na(deadwood) & deadwood == TRUE, 1, 0)
  ) %>% 
  dplyr::select(cluster, country, logging_trail, clear, grndwrk, planting, 
                anti_browsing, windthrow, deadwood) # remove unnecessary cols


## get shares : overall --------------------
df_management_shares_01 <- df_individual_management %>%
  mutate(cluster = factor(cluster)) %>%
  pivot_longer(cols = -cluster, 
               names_to = 'management_type', values_to = 'presence') %>%
  group_by(management_type) %>%
  mutate(total = n()) %>%     # Total records per management_type
  group_by(presence, management_type) %>%
  summarise(
    n = n(),
    pct = round((n / first(total)) * 100,1)  # Calculate percentage
  ) %>%
  ungroup()


# Create bar plot
p_manag_types_overall <-
  ggplot(df_management_shares_01, aes(x = management_type, y = pct, fill = factor(presence))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e")) +
  geom_text(aes(label = pct), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "white") +
  labs(
    title = "Management Types (0 vs 1)",
    x = "Management Type",
    y = "%",
    fill = ""
  ) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

# Get management shares per country ----------
## get shares : overall --------------------
df_management_shares_01_country <- df_individual_management %>%
  mutate(cluster = factor(cluster)) %>%
  pivot_longer(cols = c(-cluster, -country), 
               names_to = 'management_type', values_to = 'presence') %>%
  group_by(management_type, country) %>%
  mutate(total = n()) %>%     # Total records per management_type
  group_by(presence, management_type, country) %>%
  summarise(
    n = n(),
    pct = round((n / first(total)) * 100,1)  # Calculate percentage
  ) %>%
  ungroup()

# Define the order of countries from west to east
west_to_east_order <- c("FR", "CH", "DE", "AT",  "SI", "CZ", "SK","PL")

# Reorder the country column
df_management_shares_01_country <- df_management_shares_01_country %>%
  mutate(country = factor(country, levels = west_to_east_order))

# Plot
p_manag_types_countries <- df_management_shares_01_country %>%
  ggplot(aes(x = management_type, y = pct, fill = factor(presence))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(.~ country) +
  scale_fill_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e")) +
  geom_text(aes(label = pct), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "white") +
  labs(
    title = "Presence of Management Types (0 vs 1)",
    x = "Management Type",
    y = "%",
    fill = ""
  ) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
p_manag_types_countries

df_management_shares_01_country %>% 
  # dplyr::filter(country == 'CZ') %>% 
  ggplot(aes(x = country, y = pct, fill = factor(presence))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(.~ management_type) +
  scale_fill_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e")) +
  geom_text(aes(label = pct), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "white") +
  labs(
    title = "Presence of Management Types (0 vs 1)",
    x = "Management Type",
    y = "%",
    fill = ""
  ) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





## Get management intensity on cluster level : rescaled between 0-1 (25 is 100%, eg I divide everything by 25) ------
dat_manag_intensity_cl <- 
  dat2 %>% 
  dplyr::select(ID, cluster, manag_intensity, salvage_intensity, protection_intensity,n_plots ) %>% 
  distinct() %>% 
  group_by(cluster) %>% 
  mutate(n_plots = min(n_plots),
         manag_int_cluster  = sum(manag_intensity, na.rm =T),
         salvage_int_cluster = sum(salvage_intensity, na.rm =T),
         protect_int_cluster = sum(protection_intensity, na.rm =T)
  ) %>% 
  ungroup() %>% 
  dplyr::select(cluster, n_plots , manag_int_cluster, salvage_int_cluster,protect_int_cluster) %>% 
  distinct() %>% 
  mutate(scaled_manag_int_cluster = manag_int_cluster / (n_plots *5),
         scaled_salvage_int_cluster = salvage_int_cluster / (n_plots *3),
         scaled_protect_int_cluster = protect_int_cluster / (n_plots *2)) %>% # scale values between 0-1
  rename(management_intensity  = scaled_manag_int_cluster,
         salvage_intensity     = scaled_salvage_int_cluster,
         protection_intensity  = scaled_protect_int_cluster) %>% 
  dplyr::select(cluster, management_intensity, salvage_intensity, protection_intensity, n_plots)  #%>% 
    #View()

hist(dat_manag_intensity_cl$protection_intensity)



# regions per countries:
#     country         regions
#11   germany           11, 12, 14, 18, 19, 20, 25
#12   poland            17
#13   czech             15, 26
#14   austria           13
#15   slovakia          16
#16   slovenia          23
#17   italy             21
#18   switzerland       22
#19   france            24, 27






# "V1"
# "x_coord"
# "y_coord"
# "ID"   - for each point, gets region, group and point number
# "country"         - 9 countries   
# "region"          - 
# "group"           - cluster
# "point"           - 1-5 - check, only 5 per cluster?, max is 7
# "dist"            - disturbed? T, F - part of the Cornelius' map? 
# "clear" - cleared? T, F
# "grndwrk"        - T,F
# "logging_trail"  - T/F
# "windthrow"      - int
# "planting"       - int
# "deadwood"       - int
# "anti_browsing"
# "value"         - delete, not necessary
# "VegType" - vertical layer: regeneration, advRegeneration, Survival
# "Variable" - "n"   - count
#              "hgt" - height class
#              "dmg" - damage, only for regeneration
#              "dbh" - only for survivors
# "Species" - tree species      
# "n"        - count by species, by vertical layer

dat2 <- ungroup(dat2)
# does every point has all species?
table(dat2$ID, dat2$Species)# YES - 7 records per plot, for each species



# get summary info:
# number of points
# number of clusters

# (n_row         <- nrow(dat))
# (n_plots       <- length(unique(dat$ID)))
# (n_clusters    <- length(unique(dat$cluster)))
# (n_regions     <- length(unique(dat$region)))
# (n_countries   <- length(unique(dat$country)))


# do analyses on teh cluster evel, but considering the all vertical layers:
# are there any clusters without any regeneration?
dat2 %>% 
  group_by(cluster) %>% 
  dplyr::filter(Variable == 'n') %>% # counts, not other variables
  dplyr::filter(n==1)# %>% 



# get data for iLand: counts, species, vertical class, dbh (for Mature trees) -----
# sub_plot plot speces count VegType dbh

dat_counts <- dat2 %>% 
  dplyr::filter(Variable == 'n',
                dist == 'TRUE') %>%
  dplyr::select(ID, VegType, Species, cluster, country, n) %>% #, manag, manag_intensity,
  mutate(count = ifelse(is.na(n), 0, n)) %>% 
  dplyr::select(-n)



dat_dbh_mature <- dat2 %>%
  filter(VegType == 'Mature') %>% 
  filter(Variable == 'dbh' ) %>%  # & VegType == 'Survivor'
  dplyr::select(ID, VegType, Species,   cluster,  country, value) %>% 
  mutate(dbh = case_when(
    #VegType == 'advRegeneration' ~ 5,
    value == 4 & VegType == 'Mature' ~ 70,
    value == 3 & VegType == 'Mature' ~ 50,
    value == 2 & VegType == 'Mature' ~ 30,
    value == 1 & VegType == 'Mature' ~ 15,
    TRUE ~ NA_real_  # Default case if none of the above conditions are met
  )) #%>% 
  #dplyr::select(-value)

dat_iLand <- dat_counts %>% 
  full_join(dat_dbh_mature, by = join_by(ID, VegType, Species, cluster, country))


fwrite(dat_iLand, 'outData/sub_plots_counts_DBH.csv')


# get disturbance severity from presence/absence of mature trees:
# Remove 'Species' column and keep distinct rows

# Count unique IDs per cluster
df_subplot_counts <- dat_counts %>%
  group_by(cluster) %>%
  summarize(n_subplot = n_distinct(ID))

# summarize the count of mature trees 
df_mature_dist_severity <- dat_counts %>%
   dplyr::select(-Species) %>%
  dplyr::filter(VegType == "Mature") %>% 
  distinct() %>% 
  group_by(ID, cluster) %>% 
  summarize(count = sum(count),.groups = "drop") %>%  # Summarize count per group
  mutate(mature_presence = if_else(count > 0, 1, 0)) %>% # Add column based on count
  group_by(cluster) %>% 
  summarize(mature_presence_nplots = sum(mature_presence), .groups = "drop") %>%  # Summarize count per group
  right_join(df_subplot_counts, by = join_by(cluster)) %>% 
  mutate(mature_dist_severity = 1- mature_presence_nplots/n_subplot)

hist(df_mature_dist_severity$mature_dist_severity)


fwrite(df_mature_dist_severity, 'outData/disturb_severity_mature.csv')


#calculate stem density - convert to stems/ha ----------------------------------


# Get stem density - from vegetation matrix ------------------------------------
veg_matrix_counts <- 
  dat2 %>% 
  dplyr::filter(Variable == 'n',
                dist == 'TRUE') %>%
  dplyr::select(ID, VegType, Species, cluster, country, n) %>% #, manag, manag_intensity,
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  pivot_wider(values_from=n,
              names_from = Species) %>% 
  left_join(dat_manag_intensity_cl)


# Exclude 'ID', 'cluster', VegType columns from summarization
species_columns <- setdiff(names(veg_matrix_counts), c("ID", "cluster", 'VegType', 'country', 'management_intensity',
                                                       'salvage_intensity','protection_intensity', 'n_plots'))  #'manag', '', 

# Calculate the total and average stems per species per hectare
stem_dens_ha <-
  veg_matrix_counts %>%
  group_by(cluster, VegType,  country, management_intensity) %>% # manag, manag_intensity 
  summarize(across(all_of(species_columns), sum),
            n_plots = min(n_plots) #,
            #n_subplots = n()
            ) %>% # .groups = "drop"
    
  mutate(scaling_factor = 10000 / (4 * n_plots)) %>% # if n_subplots = 15 -ok, its accounts for 3 vertical layers
   # View()  # divide by 4 as one sub-plot = 4 m2 
    ungroup(.) %>%
  group_by(cluster, VegType, country, management_intensity) %>%  # , manag, manag_intensity
  mutate(across(all_of(species_columns), ~ .x * scaling_factor),
         total_stems_all_species = sum(across(all_of(species_columns)))) 

# 10/18/2024 corrected scalling factor
stem_dens_species_long<- 
  veg_matrix_counts %>%
  group_by(ID, cluster, VegType,  country,management_intensity ) %>% #, manag, manag_intensity
  summarize(across(all_of(species_columns), sum),
            n_plots = min(n_plots) #,
            #n_subplots = n()
            ) %>% # .groups = "drop"
  mutate(scaling_factor = 10000 / (4 * n_plots)) %>% # if n_subplots = 15 -ok, its accounts for 3 vertical layers
  #View()
    ungroup(.) %>%
  group_by(cluster, VegType, country, management_intensity) %>% # , manag, manag_intensity
  mutate(across(all_of(species_columns), ~ .x * scaling_factor),
         total_stems_all_species = sum(across(all_of(species_columns)))) %>%
  dplyr::select(-total_stems_all_species, #-n_subplots, 
                -scaling_factor ) %>%
  pivot_longer(!c(ID, cluster, VegType, country, management_intensity,n_plots ), #ID, country, management_intensity
               names_to = 'Species', 
               values_to = "stem_density") #%>%  #, manag, manag_intensity


fwrite(stem_dens_species_long, 'outData/df_stem_density_long.csv')

# calculate density per plot&species, add management intensity
stem_dens_species_long_cluster<- 
  stem_dens_ha %>% 
    ungroup(.) %>% 
   dplyr::select(-total_stems_all_species,  -scaling_factor ) %>% 
   pivot_longer(!c(cluster, VegType, country, management_intensity,n_plots  ), #ID, country, management_intensity
               names_to = 'Species', 
               values_to = "stem_density") #%>%  #, manag, manag_intensity
 

# stem density sum per vertical layer and species!
stem_dens_ha_cluster_sum <- stem_dens_ha %>% 
  group_by(cluster,  country, VegType) %>%  #
  summarize(total_stems = sum(total_stems_all_species)) #%>% 
  #dplyr::select(cluster, total_stems_all_species)
range(stem_dens_ha_cluster_sum$total_stems)

# sum stems first up across vertical groups!!
df_stems <- stem_dens_ha_cluster_sum %>% 
  group_by(cluster, country) %>% #
  summarise(sum_stems = sum(total_stems))# %>% 

range(df_stems$sum_stems)



# get stem density per species and height category - get vegetation matrixes for that - on cluster level???
# keep the NA for species as 0 - consistently across the dataset!!!
plot_density <- df_stems %>% 
    group_by(country, cluster) %>% # , Species 
  summarize(sum_n    = sum(sum_stems, na.rm = T),
            mean_n   = mean(sum_stems, na.rm = T),
            sd_n     = sd(sum_stems, na.rm = T),
            median_n = median(sum_stems, na.rm = T))


p_dens <-plot_density %>% 
  ungroup(.) %>%
  ggplot(aes(x = reorder(as.factor(country), -sum_n, FUN = median),
             y = sum_n,
             #fill = as.factor(manag),  # Use 'fill' for differentiating groups
             group = interaction(as.factor(country)))) +
  stat_summary(fun = "median", 
               geom = "bar", 
               position = position_dodge(),  # Use 'position_dodge' to place bars next to each other
               alpha = .7) +
  stat_summary(
    data = plot_density,
    mapping = aes(x = reorder(as.factor(country), -sum_n, FUN = median), 
                  y = sum_n,
                  group = interaction(as.factor(country))),
    fun.min = function(z) { quantile(z, 0.25) },
    fun.max = function(z) { quantile(z, 0.75) },
    fun = median,
    geom  = 'errorbar',
    position = position_dodge(0.9),  # Match the dodge width of the bars
    width = .2
  ) +
  xlab('') +
  ggtitle("Median stem density/cluster") + 
  theme_bw() +
  labs(fill = "", ylab = 'Stem density [/ha]') 

(p_dens)


# plot density by density plots by vertical classes


# Richness  ---------------------------------------------------------------
# Calculate species richness
# from the original table

df_richness <- 
  dat2 %>%
  dplyr::filter(Variable == 'n') %>% # counts, not other variables
  group_by(cluster,   country, Species) %>%
  summarise(sum_counts = sum(n, na.rm = T)) %>% # sum up species accros vertical groups
  #View()
    ungroup(.) %>% 
    group_by(cluster,   country) %>% 
    summarise(richness = sum((sum_counts!=0), na.rm = TRUE))
#13569
#View(df_richness)

# check summary as my plot looks a bit weird
df_richness %>% 
 # group_by(manag) %>% 
  summarise(med = median(richness, na.rm = T),
            q_25 = quantile(richness, probs = 0.25, na.rm = TRUE),
            q_75 = quantile(richness, probs = 0.75, na.rm = TRUE)
            )




# Get IVI -----------------------------------------------------------------

# spercies importance value = IVI - based on the relative density and relative basal area
# relative density: proprtion of teh species from total stems
# relative basal artea: % of species BA from all basal area
# - for BA, need to calculate the BA from the dbh values! Add for teh advanced regeneration: 2-10 cm

# get relative density:

rel_density <- 
  stem_dens_ha %>% 
  ungroup(.) %>% 
  dplyr::select( -scaling_factor, -VegType, -total_stems_all_species  ) %>%  # -n_subplots ,
  pivot_longer(piab:pist, 
               names_to = 'Species', 
               values_to = 'stem_dens') %>% 
  group_by(cluster, Species, country) %>%  #
  summarize(sum_sp_dens = sum(stem_dens, na.rm = T)) %>% # sum per species (across vertical layers)
  ungroup(.) %>%
  group_by(cluster, country) %>%  #,manag
  mutate(sum_all = sum(sum_sp_dens, na.rm = T)) %>% # sum across all species
  mutate(rel_dens = sum_sp_dens/sum_all )


# dbh table: use mean values:
# regeneration - no value, as DBH is measured at breast height
# advanced regeneraTION = 5 CM
# mature: classes: 
# 1 = 10-20cm -> 15
# 2 = 20-40cm -> 30
# 3 = 40-60cm -> 50
# 4 > 60 cm   -> 70

# add dbh based on teh vertical layer (colum Variable == 'dbh'), as advRegen and Regen do not have this value
# complete teh dbh info by vertical class
dbh_mature <- dat2 %>%
  filter(VegType == 'Mature') %>% 
  filter(Variable == 'dbh' ) %>%  # & VegType == 'Survivor'
  dplyr::select(ID, VegType, Species,   cluster,  country, value) %>% 
  mutate(dbh = case_when(
    #VegType == 'advRegeneration' ~ 5,
    value == 4 & VegType == 'Mature' ~ 70,
    value == 3 & VegType == 'Mature' ~ 50,
    value == 2 & VegType == 'Mature' ~ 30,
    value == 1 & VegType == 'Mature' ~ 15,
    TRUE ~ NA_real_  # Default case if none of the above conditions are met
  )) %>% 
  dplyr::select(-value)

dbh_mature_value <- dat2 %>%
  filter(VegType == 'Mature') %>% 
  filter(Variable == 'dbh' ) %>%  # & VegType == 'Survivor'
  dplyr::select(ID, VegType, Species,   cluster,  country, value) %>% 
  mutate(dbh = case_when(
    #VegType == 'advRegeneration' ~ 5,
    value == 4 & VegType == 'Mature' ~ 70,
    value == 3 & VegType == 'Mature' ~ 50,
    value == 2 & VegType == 'Mature' ~ 30,
    value == 1 & VegType == 'Mature' ~ 15,
    TRUE ~ NA_real_  # Default case if none of the above conditions are met
  )) #%>% 
  #dplyr::select(-value)

  
dbh_advanced <- 
  dat2 %>%
  filter(VegType == 'Juveniles') %>% 
  filter(Variable == 'n') %>% 
  dplyr::select(ID, VegType, Species,   cluster,  country, value) %>% 
  mutate(dbh = case_when(!is.na(value) & VegType == 'Juveniles' ~ 5,
    TRUE ~ NA_real_  # Default case if none of the above conditions are met
  )) %>% 
  dplyr::select(-value)

# account for the regeneration as basal area: use value of 1 mm
dbh_regen <- 
  dat2 %>%
  filter(VegType == 'Saplings') %>% 
  filter(Variable == 'n') %>% 
  dplyr::select(ID, VegType, Species,   cluster,  country, value) %>% 
  #mutate(dbh = NA)  %>% 
  mutate(dbh = case_when(!is.na(value) & VegType == 'Saplings' ~ 0.01,  # 0.1 mm
                         TRUE ~ NA_real_  # Default case if none of the above conditions are met
  )) %>% 
  dplyr::select(-value)


# merge 3 vertical classes
dbh_merge_vert <- bind_rows(dbh_mature, dbh_advanced, dbh_regen)

#dbh_merge_vert_cluster <- dbh_merge_vert %>% 
#  group_by(cluster)

df_out_density_dbh <- stem_dens_species_long %>% 
  left_join(dbh_merge_vert, 
            by = c('ID', 'cluster', 'VegType', 'country', 'Species')) %>%  # , 'manag'
  group_by(cluster, VegType, country, management_intensity, Species, dbh) %>% 
  summarise(stem_density = sum(stem_density, na.rm = T)) %>% 
  ungroup()

table(df_out_density_dbh$Species, df_out_density_dbh$VegType)
# export dbh for iLand? maybe it is later needed?
fwrite(df_out_density_dbh, 'outData/veg_density_DBH.csv')

# add DBH info to the stem density one to calculate Basal area
df_BA <- 
  stem_dens_species_long %>% 
  left_join(dbh_merge_vert, by = c('ID', 'cluster', 'VegType', 'country', 'Species')) %>%  # , 'manag'
  mutate(r = dbh/2,
         BA = pi*r^2,
         BA_stems = BA*stem_density) %>% 
  ungroup(.) %>% 
  dplyr::select(-VegType) %>%
  group_by(cluster,Species,   country) %>% 
    summarise(BA_sp = sum(BA_stems, na.rm = T)) %>% 
    ungroup(.) %>% 
    group_by(cluster) %>% 
    mutate(BA_sum = sum(BA_sp, na.rm = T),
           rel_BA = BA_sp/BA_sum)
  
#View(df_BA)

# calculate species importance value: based on relative density and relative basal area per species
df_IVI <- 
  rel_density %>% 
  full_join(df_BA, by = join_by(cluster, Species, country)) %>% # manag 
  mutate(rIVI = ( rel_dens +rel_BA)/2) %>%  # relative IVI
  tidyr::replace_na(., list(rIVI = 0, rel_BA   = 0)) 



# Vertical layers ---------------------------------------------------------
df_vert <- 
  stem_dens_species_long %>% 
  dplyr::filter(stem_density > 0) %>% 
  ungroup(.) %>% 
  dplyr::select(cluster, country,   VegType) %>%
  group_by(cluster) %>% # , manag
  distinct(.) %>% 
  dplyr::summarize(n_layers = n()) #%>% 

# amke sure that all clusters are present, add 0s oif there is no stems found
df_vert <- merge(df_stems, df_vert, by = c("cluster"), all.x = TRUE) # 
df_vert$n_layers[is.na(df_vert$n_layers)] <- 0




# export data: -----------------------------------------------------------------
save(df_individual_management,   # individual management types (0,1) per subplot 
     dat_manag_intensity_cl, # get the management intensity value per cluster
     df_IVI,  # species importance value
     df_richness, 
     df_vert,   # number of vertical layers
     df_stems,  # numbers of stems for all vert layers
     stem_dens_ha_cluster_sum,   
     stem_dens_species_long,  # stem density per ID& species
     stem_dens_species_long_cluster, # summaziesd per cluster, species and vert class
     veg_matrix_counts,    # stem counts (not stem density)
     file="outData/veg.Rdata")


