

# extract variables 
# density
# tree species
# IVI
# how to summarize data? average oevr clusters???
# smart way to keep empty clusters?
# make a major df


# analyse structure and composition
# recalculate to hectares


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

dat <- fread('rawData/working_directory/rapid_assessment_mdf.csv')

head(dat)

summary(dat)
names(dat)

# data

dim(dat)


# correct management based on Christians scrips
dat2 <- dat %>% 
  mutate(
    logging_trail = ifelse(!is.na(logging_trail) & logging_trail == "true", 1, 0),
    clear = ifelse(!is.na(clear) & clear == "true", 1, 0),
    grndwrk = ifelse(!is.na(grndwrk) & grndwrk == "true", 1, 0),
    planting = ifelse(!is.na(planting) & planting != 2, 1, 0),
    anti_browsing = ifelse(!is.na(anti_browsing) & anti_browsing != 2, 1, 0),
    management = case_when(
      logging_trail == 1 | clear == 1 | grndwrk == 1 | planting == 1 | anti_browsing == 1 ~ 1,
      TRUE ~ 0  # Default case
    )
  )




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
#              "hgt" - height?
#              "dmg" - damage, only for regeneration
#              "dbh" - only for survivors
# "Species" - tree species      
# "n"        - count by species, by vertical layer


# Classes ‘data.table’ and 'data.frame':	1264830 obs. of  21 variables:
#   $ V1           : int  1 2 3 4 5 6 7 8 9 10 ...
# $ x_coord      : int  1245084 1245113 1245084 1245085 1245050 1245152 1245154 1245121 1245186 1245153 ...
# $ y_coord      : int  6744529 6744529 6744498 6744562 6744527 6744883 6744852 6744883 6744881 6744915 ...
# $ ID           : chr  "11_11_101_1" "11_11_101_2" "11_11_101_3" "11_11_101_4" ...
# $ country      : int  11 11 11 11 11 11 11 11 11 11 ...
# $ region       : int  11 11 11 11 11 11 11 11 11 11 ...
# $ group        : int  101 101 101 101 101 102 102 102 102 102 ...
# $ point        : int  1 2 3 4 5 1 2 3 4 5 ...
# $ dist         : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
# $ clear        : logi  TRUE NA TRUE TRUE TRUE TRUE ...
# $ grndwrk      : logi  TRUE NA TRUE TRUE TRUE NA ...
# $ logging_trail: logi  NA TRUE NA NA NA NA ...
# $ windthrow    : int  3 3 3 3 3 3 3 3 3 3 ...
# $ planting     : int  3 3 3 3 3 2 3 2 3 2 ...
# $ deadwood     : int  NA NA NA NA NA NA NA NA NA NA ...
# $ anti_browsing: int  3 3 3 3 3 3 3 3 3 3 ...
# $ value        : chr  NA NA NA NA ...
# $ VegType      : chr  "Regeneration" "Regeneration" "Regeneration" "Regeneration" ...
# $ Variable     : chr  "n" "n" "n" "n" ...
# $ Species      : chr  "piab" "piab" "piab" "piab" ...

dat %>% 
  dplyr::filter(country == 11 & region == 11& group == 102 ) %>% 
  dplyr::filter(!is.na(n)) %>% 
  dplyr::select(c(ID, country, group, point, region, value, VegType, Variable, Species, n)) %>%
  arrange(point)
  


# check if there are any clusters without any vegetation on them:
dat %>% 
group_by(country, group, point) %>% 
dplyr::filter(is.na(n)) %>% 
  dplyr::select(c(ID, country, group, point, value, VegType, Variable, Species, n)) %>%
  arrange(point)


# get density per plot and species - create vegetation matrixes

data <- data.frame(
  plot = c('A', 'A', 'A', 'B', 'B'),
  sub_plot = c(1, 1, 1, 1, 1),
  type = c(NA, 'manag', NA, NA, 'unmanag')
)

# Replace NA by the non-NA value in the group
#data <- 
  data %>%
  group_by(plot, sub_plot) %>%
  mutate(type = ifelse(is.na(type), first(type[!is.na(type)]), type)) %>%
  ungroup()



# try t replace NA values to cefine the managed sites: needs some more work!
  # eg make sure that each ID has at least one management indication(i.e. planted, ground work, fencing...)
dat <- dat %>% 
  mutate(cluster = paste(region, group, sep = '_'), # create unique ID per cluster
         area = 4)  %>%                       # add recording area: 4 m2 to calculae the density for hectar
  group_by(ID) %>% 
 # mutate(clear = ifelse(is.na(clear), 
#                       max(clear, na.rm = TRUE),  # Use max for logical
#                       clear)) %>%
  ungroup(.)
  
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








# get stem density per species and height category - get vegetation matrixes for that - on cluster level???
# keep the NA for species as 0 - consistently across the dataset!!!
df_density_as0 <- dat %>% 
  dplyr::filter(Variable == 'n') %>% 
  mutate(n = ifelse(is.na(n), 0, n)) %>% 
  group_by(country, cluster, point, VegType, Species ) %>% 
  summarize(sum_n = sum(n, na.rm = T))

# plot: what is the density of regeneration per plot & country?
df_density_as0 %>% 
  ungroup(.) %>% 
  filter(VegType == 'Regeneration') %>% 
  dplyr::select(-Species) %>% 
  filter(sum_n !=0) %>% 
  ggplot(aes(x = as.factor(country),
             y = sum_n,
             group = country)) +
  geom_boxplot()
  
  # get average per cluster?
  


