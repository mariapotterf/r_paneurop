

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
# 


# ask Christinan - included empty plots?
# what is 'cluster?  = group?
# country = '

library(data.table)
library(dplyr)

dat <- fread('rawData/working_directory/rapid_assessment_mdf.csv')

head(dat)

summary(dat)
names(dat)

# data

dim(dat)

# "V1"
# "x_coord"
# "y_coord"
# "ID"   - for each point, gets region, group and point number
# "country"         - 9 countries   
# "region"          - 
# "group"           - cluster
# "point"           - 1-5 - check, only 5 per cluster?, max is 7
# "dist"            - disturbed? T, F - part of the Cerneliuses map? 
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
dplyr::filter(is.na(n)) #%>% 
  dplyr::select(c(ID, country, group, point, value, VegType, Variable, Species, n)) %>%
  arrange(point)


# get 

unique(dat$Species)  # 

#[1] "piab" "pisy" "lade" "abal" "psme" "taba" "fasy" "quro" "acca" "acpl" "acps" "algl" "alin" "alvi"
#[15] "potr" "posp" "besp" "frex" "tisp" "prav" "soau" "soto" "soar" "casa" "aehi" "cabe" "ulsp" "rops"
#[29] "saca" "juni" "jure" "qusp" "sasp" "aial" "osca" "fror" "ilaq" "pist"

# get density per plot and species - create vegetation matrixes

dat <- dat %>% 
  mutate(cluster = paste(region, group, sep = '_'), # create unique ID per cluster
         area = 4)                        # add recording area: 4 m2 to calculae the density for hectar

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


# get stem density per species and height category
df_density <- dat %>% 
  group_by()
  


