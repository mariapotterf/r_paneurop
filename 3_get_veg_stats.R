

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
# IVI - species importance value, need to get BA
# # include management


# ask Christinan - included empty plots? YES
# what is 'cluster?  = region +group?
# country = '
# management = fill in missing values correctly: based on the first estimation per plot (without completed species)

library(data.table)
library(dplyr)
library(ggplot2)
library(terra)
library(tidyr)


dat <- fread('rawData/working_directory/rapid_assessment_mdf.csv')

#table(dat$dist)

# create management based on Christians scrips
dat<- dat %>% 
  filter(dist == TRUE) %>% # remove the plot if not disturbed
  mutate(cluster = paste(region, group, sep = '_')) %>% 
  mutate(
    logging_trail = ifelse(!is.na(logging_trail) & logging_trail == "true", 1, 0),
    clear         = ifelse(!is.na(clear) & clear == "true", 1, 0),
    grndwrk       = ifelse(!is.na(grndwrk) & grndwrk == "true", 1, 0),
    planting      = ifelse(!is.na(planting) & planting != 2, 1, 0),
    anti_browsing = ifelse(!is.na(anti_browsing) & anti_browsing != 2, 1, 0),
    manag    = case_when(
      logging_trail == 1 | clear == 1 | grndwrk == 1 | planting == 1 | anti_browsing == 1 ~ 'Managed',
      TRUE ~ 'Unmanaged'  # Default case
    )
  )

dat <- dat %>% 
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




# Create master df with empty plots - eg no trees found on them
df_master <- 
  dat %>% 
  dplyr::select(country,region, group, cluster, point, ID) %>%  #  region, group, 
    unique() %>%  # remove duplicated rows
  group_by(country, region, group, cluster) %>%  # region, group,
  summarize(n_plots = n())

# filter only clusters that have 4 points

# keep clusters with 4, as some disturbace plots were very small
# remove if there is less records
dat <- dat %>% 
  right_join(df_master, by = join_by(country, region, group, cluster)) %>%
  filter(n_plots > 3)
  





#View(df_master)
head(dat)

summary(dat)
names(dat)

# data

dim(dat)

#11 germany
#12 poland
#13 czech
#14 austria
#15 slovakia
#16 slovenia
#17 italy
#18 switzerland
#19 france


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
#19   france            24






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
#"anti_browsing"
# "value"         - delete, not necessary
# "VegType" - vertical layer: regeneration, advRegeneration, Survival
# "Variable" - "n"   - count
#              "hgt" - height class
#              "dmg" - damage, only for regeneration
#              "dbh" - only for survivors
# "Species" - tree species      
# "n"        - count by species, by vertical layer

dat %>% 
  dplyr::filter(country == 'DE' & region == 11& group == 102 ) %>% 
  dplyr::filter(!is.na(n)) %>% 
  dplyr::select(c(ID, country, group, point, region, value, VegType, Variable, Species, n)) %>%
  arrange(point)
  


# check if there are any clusters without any vegetation on them:
dat %>% 
group_by(country, group, point) %>% 
dplyr::filter(is.na(n)) %>% 
  dplyr::select(c(ID, country, group, point, cluster, manag, VegType, Variable, Species, n)) %>%
  arrange(point)


# does every point has all species?
table(dat$ID, dat$Species)# YES - 7 records per plot, for each species

# get summary info:
# number of points
# number of clusters

(n_row         <- nrow(dat))
(n_plots       <- length(unique(dat$ID)))
(n_clusters    <- length(unique(dat$cluster)))
(n_regions     <- length(unique(dat$region)))
(n_countries   <- length(unique(dat$country)))


# do analyses on teh cluster evel, but considering the all vertical layers:
# are there any clusters without any regeneration?
dat %>% 
  group_by(cluster) %>% 
  dplyr::filter(Variable == 'n') %>% # counts, not other variables
  dplyr::filter(n==1)# %>% 

# seems that all clusters have at least one stem left


# Dummy example: calculate stem density - convert to stems/ha ----------------------------------
# accont for missing species

# calculate from vegetation matrix 
dd <- data.frame(cluster = c(rep('a',5),rep('b',3)),
                 ID = c(1:5, 1:3),
                 sp1 = c(0, 1,1,1,1, 0,0,0),
                 sp2 = c(0, 0,5,0,0,0,0,1),
                 sp3 = c(1, 0,0,0,4,0,0,16),
                 sp4 = c(0, 0,0,0,0,0,0,0))

(dd)

# Exclude 'ID' and 'cluster' columns from summarization
species_columns <- setdiff(names(dd), c("ID", "cluster"))

# Calculate the total and average stems per species per hectare
dd_summary <- 
  dd %>%
  group_by(cluster) %>%
  summarize(across(all_of(species_columns), sum),
            rows_per_cluster = n()) %>% # .groups = "drop"
  mutate(scaling_factor = 10000 / (4 * rows_per_cluster)) %>%
    group_by(cluster) %>% 
    mutate(across(all_of(species_columns), ~ .x * scaling_factor),
         total_stems_all_species = sum(across(all_of(species_columns))))

dd_summary


# Get stem density - from vegetation matrix ------------------------------------
veg_matrix_counts <- 
  dat %>% 
  dplyr::filter(Variable == 'n',
                dist == 'TRUE') %>%
  dplyr::select(ID, VegType, Species, cluster, manag, country, n) %>% 
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  pivot_wider(values_from=n,names_from = Species)


# Exclude 'ID', 'cluster', VegType columns from summarization
species_columns <- setdiff(names(veg_matrix_counts), c("ID", "cluster", 'VegType', 'manag', 'country'))

dat %>% 
  filter(cluster == '15_110') %>%
  View()
  #distinct(point)

length(unique(veg_matrix_counts$cluster))
length(unique(stem_dens_ha_cluster_sum$cluster))

# Calculate the total and average stems per species per hectare
stem_dens_ha <- veg_matrix_counts %>%
  group_by(cluster, VegType,  country, manag ) %>%
  summarize(across(all_of(species_columns), sum),
            rows_per_cluster = n()) %>% # .groups = "drop"
  mutate(scaling_factor = 10000 / (4 * rows_per_cluster)) %>% # if rows_per_cluster = 15 -ok, its accounts for 3 vertical layers
  ungroup(.) %>% 
  group_by(cluster, VegType, country, manag) %>% 
  mutate(across(all_of(species_columns), ~ .x * scaling_factor),
         total_stems_all_species = sum(across(all_of(species_columns))))

# stem density table
stem_dens_ha_cluster <- stem_dens_ha %>% 
  dplyr::select(cluster, country, manag, VegType, total_stems_all_species)

stem_dens_ha_cluster_sum <- stem_dens_ha %>% 
  group_by(cluster,  country, manag, ) %>% 
  summarize(sum_stems = sum(total_stems_all_species)) #%>% 
  #dplyr::select(cluster, total_stems_all_species)


View(stem_dens_ha)
View(stem_dens_ha_cluster)

# cluster 21_143 has 21 of teh records!! - keep, 7 points in cluster
dat %>% 
  filter(cluster == '21_143') %>% 
  dplyr::select(ID, VegType, dist, n) %>% 
  distinct() %>% 
  View()


# get stem density per species and height category - get vegetation matrixes for that - on cluster level???
# keep the NA for species as 0 - consistently across the dataset!!!
plot_density <- stem_dens_ha_cluster %>% 
    group_by(country, manag, cluster) %>% # , Species 
  summarize(sum_n    = sum(total_stems_all_species, na.rm = T),
            mean_n   = mean(total_stems_all_species, na.rm = T),
            sd_n     = sd(total_stems_all_species, na.rm = T),
            median_n = median(total_stems_all_species, na.rm = T))


df_density_as0 %>% 
  ungroup(.) %>%
  ggplot(aes(x = reorder(as.factor(country), -sum_n, FUN = median),
             y = sum_n,
             fill = as.factor(manag),  # Use 'fill' for differentiating groups
             group = interaction(as.factor(country), as.factor(manag)))) +
  stat_summary(fun = "median", 
               geom = "bar", 
               position = position_dodge(),  # Use 'position_dodge' to place bars next to each other
               alpha = .7) +
  stat_summary(
    data = df_density_as0,
    mapping = aes(x = reorder(as.factor(country), -sum_n, FUN = median), 
                  y = sum_n,
                  group = interaction(as.factor(country), as.factor(manag))),
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
  labs(fill = "") 


# need to accounto for zeros here!!!!------------------------------------------ 
mean(c(0,0,0,0,5,6,1))
mean(c(5,6,1))
median(c(0,0,0,0,5,6,1))
median(c(5,6,1))

# check how not specifically calculating for 0 affects cluster results

# 

# Richness  ---------------------------------------------------------------
df_richness <- dat %>% 
  dplyr::filter(Variable == 'n') %>% # counts, not other variables
  # filter only plots with some density
  dplyr::filter(n>0) %>% 
  group_by(country, cluster, manag) %>% 
  # count number of individual Species per group
  distinct(Species) %>% 
  count() %>%
  rename(richness = n)



# check summary as my plot looks a bit weird
df_richness %>% 
  group_by(country, manag) %>% 
  summarise(med = median(richness, na.rm = T),
            q_25 = quantile(richness, probs = 0.25, na.rm = TRUE),
            q_75 = quantile(richness, probs = 0.75, na.rm = TRUE)
            )



df_richness %>% 
  ungroup(.) %>%
  ggplot(aes(x = reorder(as.factor(country), -richness, FUN = median),
             y = richness,
             fill = as.factor(manag),  # Use 'fill' for differentiating groups
             group = interaction(as.factor(country), as.factor(manag)))) +
  stat_summary(fun = "median",
               geom = "bar",
               position = position_dodge(),  # Use 'position_dodge' to place bars next to each other
               alpha = .7) +
  stat_summary(
    data = df_richness,
    mapping = aes(x = reorder(as.factor(country), -richness, FUN = median), 
                  y = richness,
                  group = interaction(as.factor(country), as.factor(manag))),
    fun.min = function(z) { quantile(z, 0.25) },
    fun.max = function(z) { quantile(z, 0.75) },
    fun = median,
    geom  = 'errorbar',
    position = position_dodge(0.9),  # Match the dodge width of the bars
    width = .2
  ) +
  xlab('') +
  ggtitle("Richness (per cluster)") + 
  theme_bw() +
  labs(fill = "") 



# Merge structure & composition into single space -------------------------

df_org_space <- df_density_as0 %>% 
  full_join(df_richness, by = join_by(country, manag, cluster))


# get summary for the scatter plot, one point is one country & management
df_summary <- df_org_space %>%
  group_by(country, manag) %>%
  summarize(
    rich_mean = mean(richness, na.rm = TRUE),
    rich_sd = sd(richness, na.rm = TRUE),
    rich_min = rich_mean -rich_sd, #quantile(richness, 0.25, na.rm = TRUE),
    rich_max = rich_mean +rich_sd,# quantile(richness, 0.75, na.rm = TRUE),
    dens_mean = mean(sum_n, na.rm = TRUE),
    dens_sd = sd(sum_n, na.rm = TRUE),
    dens_min = dens_mean - dens_sd, #quantile(sum_n, 0.25, na.rm = TRUE),
    dens_max = dens_mean + dens_sd, #quantile(sum_n, 0.75, na.rm = TRUE),
    .groups = 'drop'
  )




df_org_space %>% 
  ggplot(aes(x = richness,
         y = sum_n,
         color = factor(manag)))+ 
  geom_point() +
  facet_wrap(.~country)




library(RColorBrewer)
my_colors <- brewer.pal(5, "Greens")

# Add 'white' at the beginning
my_colors <- c("white", my_colors)




ggplot(df_org_space) +
  stat_density_2d(aes(x = richness, y = sum_n, fill = after_stat(-level)), 
                  geom = "polygon",
                  alpha = 0.6) +
  geom_point(data = df_summary, 
             aes(x = rich_mean, y = dens_mean, shape = factor(manag), color = factor(manag)), 
             size = 1.5) +
  scale_color_manual(values = c("red", 'black')) +  # Shapes for two management categories (e.g., circle and triangle)
  
  scale_shape_manual(values = c(16, 17)) +  # Shapes for two management categories (e.g., circle and triangle)
  geom_errorbar(data = df_summary, 
                aes(x = rich_mean, 
                    ymin = dens_min, 
                    ymax = dens_max, 
                    color = factor(manag)),
                width = 0.3, 
                linewidth = 0.2) +
  geom_errorbarh(data = df_summary, 
                 aes(y = dens_mean, 
                     xmin = rich_min, 
                     xmax = rich_max, 
                     color = factor(manag)),
                 height = 2, 
                 width = 2,
                 linewidth = 0.2) +
  ylab('Mean stem density') +
  xlab(bquote('Mean richness')) +
  theme_minimal() +
  facet_wrap(.~country, scales = 'free') + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        aspect.ratio = 1,
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_blank()) #+
  #scale_color_brewer(palette = "Set1")  




ggplot(df_org_space) +
 # stat_density_2d(aes(x = richness, y = sum_n, fill = after_stat(-level)), 
  #                geom = "polygon",
   #               alpha = 0.6) +
  geom_jitter(#data = df_summary, 
             aes(x = richness, y = sum_n, shape = factor(manag), color = factor(manag)), 
             size = 1,
             alpha = 0.5) +
  scale_color_manual(values = c("red", 'black')) +  
  scale_shape_manual(values = c(16, 17)) +  # Shapes for two management categories (e.g., circle and triangle)
  ylab('Mean stem density') +
  xlab(bquote('Mean richness')) +
  theme_minimal() +
  facet_grid(.~manag) +
  theme(legend.position = 'none',
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        aspect.ratio = 1,
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_blank()) #+
#scale_color_brewer(palette = "Set1")  



