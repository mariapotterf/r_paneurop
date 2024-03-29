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

library(ggpubr)

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

library(AER)  # for hurdle model


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

# get only samplings
df_stems_sapl <- stem_dens_ha_cluster_sum %>% 
  filter(VegType == 'Regeneration') %>% 
  group_by(cluster, manag) %>% 
  summarise(sum_stems = sum(total_stems))# %>% 

# get saplings and juveniles
df_stems_both <- stem_dens_ha_cluster_sum %>% 
  filter(VegType != 'Survivor') %>% 
  group_by(cluster, manag) %>% 
  summarise(sum_stems = sum(total_stems))# %>% 

table(stem_dens_species_long$VegType,stem_dens_species_long$Species)


# aggregate by cluster: sum up teh species, by height category
stem_dens_species_cluster_long_share <- stem_dens_species_long %>% 
  ungroup(.) %>% 
  dplyr::select(-manag) %>% 
  group_by(cluster, Species, VegType) %>%
    summarise(stem_density = sum(stem_density))
    

fwrite(stem_dens_species_cluster_long_share, 'vegData_cluster_densities.csv')

disturbance_chars <- disturbance_chars %>%
  mutate(country = as.character(country))

# merge all predictors together, ID level
df_predictors <- 
  climate %>% 
  left_join(dplyr::select(distance_to_edge, c(-country, -dist_ID ))) %>% 
  left_join(dplyr::select(disturbance_chars, c(-region, -country ))) %>% 
  left_join(soil) %>%
  left_join(dplyr::select(terrain, c(-country, -region))) %>% 
  mutate(cluster = str_sub(ID, 4, -3)) 

#fwrite(df_predictors, 'outData/all_predictors.csv')

# select the dominant species per cluster
df_IVI_max <- df_IVI %>% 
  dplyr::select(cluster, Species, manag, country,rIVI) %>% 
  ungroup(.) %>% 
  group_by(cluster, country, manag) %>% 
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
  dplyr::select(cluster, manag) %>% 
  group_by(cluster, manag) %>% 
  distinct(.) %>%
  ungroup(.) %>% 
  group_by(cluster) %>% 
  summarise(type = case_when(
    all(manag == 'Unmanaged') ~ 'Unmanaged',
    all(manag == 'Managed') ~ 'Managed',
    TRUE ~ 'Mixed'
  )) %>%
  ungroup(.) %>% 
  dplyr::select(-cluster) %>% 
  table()


(manag_types)  # 849

#MAN MIX UNM 
#616 149  84 

# remove 'MIXED' - as they likely have only less plots

cluster_id_keep <- stem_dens_species_long %>%
  ungroup(.) %>% 
  dplyr::select(cluster, manag) %>% 
  group_by(cluster, manag) %>% 
  distinct(.) %>%
  ungroup(.) %>% 
  group_by(cluster) %>% 
  summarise(type = case_when(
    all(manag == 'Unmanaged') ~ 'Unmanaged',
    all(manag == 'Managed') ~ 'Managed',
    TRUE ~ 'Mixed'
  )) %>%
  dplyr::filter(type != 'Mixed') %>% 
  pull(cluster)
  #ungroup(.) %>% 
  #dplyr::select(-cluster) %>% 
  #table()




# 2. species composition per stems: saplings vs juveniles?? --------------
# how many plots with Mature trees?
# 
dom_species <- stem_dens_species_long %>% 
  dplyr::filter(cluster %in% cluster_id_keep) %>% 
  filter(VegType != 'Survivor' ) %>%  # Survivors: on 191
  group_by(Species, VegType, manag ) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup(.) %>% 
  group_by(VegType,manag) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share = sum_stems/sum_vegType*100) #%>% 

# dominant species by country?
dom_species_country <- 
  stem_dens_species_long %>% 
  dplyr::filter(cluster %in% cluster_id_keep) %>% 
  filter(VegType != 'Survivor' ) %>%  # Survivors: on 191
  group_by(Species, VegType, manag, country ) %>% 
  summarize(sum_stems = sum(stem_density, na.rm = T)) %>% 
  ungroup(.) %>% 
  group_by(VegType,manag, country) %>% 
  mutate(sum_vegType = sum(sum_stems),
         share = sum_stems/sum_vegType*100)



# SElect the 5 species with teh highest share per management, vertical structure
top_species_per_country <- dom_species_country %>%
  ungroup(.) %>% 
  group_by(country, VegType, manag) %>%
  arrange(desc(share)) %>%  # Arrange in descending order of share within each group
  slice_max(order_by = share, n = 5)  # Select top 5 species with max share in each group

View(top_species_per_country)



# how many species in total?
distinct_species_per_manag <- dom_species %>%
  filter(sum_stems > 0) %>%
  group_by(manag) %>%
  summarise(distinct_species_count = n_distinct(Species))

distinct_species_per_manag

# SElect the 5 species with teh highest share per management, vertical structure
top_species_per_group <- dom_species %>%
  group_by(VegType, manag) %>%
  arrange(desc(share)) %>%  # Arrange in descending order of share within each group
  slice_max(order_by = share, n = 5)  # Select top 5 species with max share in each group

(top_species_per_group)

# see all species
top_species_per_group %>% 
  #filter(
  #dplyr::filter(share > 1) %>% 
  ggplot(aes(x = manag,
             y = share,
             group = VegType,
             fill = Species)) +
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1)) +
  facet_grid(.~VegType)
  
  #geom_col('')


# Reorder Species within each group based on share
top_species_per_group <- top_species_per_group %>%
  arrange(VegType, manag, desc(share)) %>%
  mutate(Species = factor(Species, levels = unique(Species)),
         VegType = factor(VegType, 
                          levels = c("Regeneration", "advRegeneration"),
                          labels = c("Saplings", "Juveniles")))


# Create a bar plot with ordered categories
ggplot(top_species_per_group, aes(x = manag, y = share, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(share), y = share), 
            position = position_dodge(width = 0.9), 
            color = "black", size = 3, vjust = -0.5) + # Adjust vjust for label position
  facet_grid(. ~ VegType) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) +
  labs(y = "Share [%]", x = "", fill = "Species") +
  theme_classic()


ggplot(top_species_per_group, aes(x = manag, y = share, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(#aes(label = round(share,1)), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3) +
  facet_grid(. ~ VegType) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1)) +
  labs(y = "Share [%]", x = "Management Type", fill = "Species")


# all species by aboundance
dom_species %>% 
  dplyr::filter(share > 1) %>% 
  ggplot(aes(x = manag,
             y = share,
             fill = Species)) +
  geom_bar(position = "fill", stat = "identity")+
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
  facet_grid(VegType~manag)
#geom_col('')


# 3. how many clusters/plots has less then 2000 regen/ha? ----------------------------------
# which country is the most affected? (where they ccur most often?)



df_stock_density <- 
  df_stems %>% 
  dplyr::filter(cluster %in% cluster_id_keep) %>% 
  mutate(dens_category = case_when(sum_stems == 0 ~ '0',
                                   sum_stems >= 1   & sum_stems < 500 ~ '1-500',
                                   sum_stems >= 500 & sum_stems < 1000 ~ '500-1000',
                                  # sum_stems >= 1000 & sum_stems < 1500 ~ '1000-1500',
                                   sum_stems >= 1000 & sum_stems < 2000 ~ '1000-2000',
                                   sum_stems >= 2000 ~ '>2000'
                                   )) %>%
  mutate(dens_category_simpl = case_when(sum_stems < 1000 ~ '<1000',
                                         sum_stems >= 1000 &sum_stems < 2000 ~ '1000-2000',
                                         sum_stems >= 2000 ~ '>2000'
  )) %>%
  mutate(dens_category = factor(dens_category, levels = 
                                  c('0', '1-500','500-1000','1000-2000', '>2000' ))) %>% 
  mutate(dens_category_simpl = factor(dens_category_simpl, levels = 
                                c('<1000','1000-2000', '>2000' ))) 
str(df_stock_density)

df_stock_density %>% 
  ggplot(aes(x = manag,
             fill = dens_category)) +
  geom_bar(position = 'fill') +
  scale_fill_manual(values = c("0" = "red", 
                               "500-1000" = "lightgreen", 
                               "1000-2000" = "green", 
                               ">2000" = "darkgreen",
                               "1-500" = "forestgreen")) +
  ylab('Share [%]') + 
  theme_classic() +
  labs(fill = "Density class") +
  xlab('')

unique(df_stock_density$dens_category)

#### Pearson chi squared test --------------------------
# Create a contingency table
contingency_table <- table(df_stock_density$manag, df_stock_density$dens_category)

# Perform the chi-squared test
chi_squared_test <- chisq.test(contingency_table)

# View the results of the chi-squared test
print(chi_squared_test)

# invalid chi suqre test - 0 categories, very little counts (less then 5)
# Combine sparse categories
df_stock_density$combined_category <- ifelse(df_stock_density$dens_category %in% c("0", "1-500"), 
                                             "0-1-500", 
                                             df_stock_density$dens_category)

# Create a new contingency table with combined categories
new_contingency_table <- table(df_stock_density$manag, df_stock_density$combined_category)

# Perform the chi-squared test on the new table
new_chi_squared_test <- chisq.test(new_contingency_table)

# View the results
print(new_chi_squared_test)



# Calculate the percentages
df_stock_density_simpl <- df_stock_density %>%
  ungroup(.) %>% 
  count(manag, dens_category) %>%
  group_by(manag) %>%
  mutate(freq = n / sum(n),
         label = scales::percent(freq, accuracy = 0.1))


df_stock_density_simpl %>% 
  ggplot( aes(x = manag, y = freq, fill = dens_category)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = label, y = freq/2), position = position_fill(vjust = 0.5), size = 4) +
  ylab('Share [%]') +
  theme_classic() +
  labs(fill = "Density class")



  

# Compare categories occurences by chi square


# How many vertical layers I have where? --------------------------------------

# Get again individual layers:
df_vert_full <- 
  stem_dens_species_long %>% 
  filter(cluster %in% cluster_id_keep) %>%  
  dplyr::filter(stem_density>0) %>% 
  ungroup(.) %>% 
  dplyr::select(cluster, country, manag, VegType) %>%
  group_by(cluster, country, manag) %>%
  distinct(.) %>% 
  mutate(n_layers = n()) #%>% 



# Calculate frequencies of individual layers and combinations
layer_frequencies <- 
  df_vert_full %>%
  group_by(manag, cluster) %>% 
  right_join(df_stems) %>% # to account for teh empty ones
    dplyr::filter(cluster %in% cluster_id_keep) %>% 
   # View()
    mutate(VegType = case_when(
      is.na(VegType) ~ 'NO',   # Replace NA with 'NO'
      TRUE ~ VegType           # Keep existing values for non-NA cases
    )) %>% 
  summarise(layers = paste(sort(unique(VegType)), collapse = "-")) %>%
  ungroup() %>%
  count(manag, layers)


# View the result
print(layer_frequencies)

ggplot(layer_frequencies, aes(x = manag, y = n, fill = layers)) +
  geom_bar(stat = "identity", position = "fill") +
  #geom_text(aes(label = paste0(round(scales::percent, 1), "%"), 
  #              y = cumsum(n) - n/2), 
  #          position = position_fill(), 
  #          color = "white", 
  #          size = 3) +
  ylab("Percentage (%)") +
  theme_classic() +
  labs(fill = "Layers")


# 4. stems vs weather - mean 2019-2022, or SPEI ---------------------------------
# average predictors by cluster
df_predictors_clust <- df_predictors %>% 
  dplyr::filter(cluster %in% cluster_id_keep) %>% 
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


# get cluster number
spei <- spei %>% 
  dplyr::mutate(year  = lubridate::year(date),
                month = lubridate::month(date)) %>%  
  dplyr::select(-c(date, -scale)) %>% 
  mutate(cluster = str_sub(ID, 4, -3))#%>%

spei_sub <- spei %>% 
    ungroup(.) %>% 
    dplyr::filter(year %in% 2018:2020 & month %in% 4:9) %>% # select just vegetation season
  dplyr::filter_all(all_vars(!is.infinite(.))) %>% # remove all infinite values
  #View()
    group_by(cluster) %>%
    summarise(spei = mean(spei, na.rm = T)) #%>% 
    
    
# join to stems counts: filter, as from veg data I have already excluded Italy and teh uncomplete clusters
df_fin <- #df_stems_both %>%
  df_stems %>% 
  left_join(df_predictors_clust, by = join_by(cluster)) %>% 
  left_join(spei_sub, by = join_by(cluster)) %>% 
  mutate(manag = as.factor(manag)) %>% 
  mutate(sum_stems = round(sum_stems)) %>% 
  dplyr::filter(cluster %in% cluster_id_keep) %>% 
  na.omit() %>% 
  as.data.frame()


length(unique(df_fin$cluster))

# clusters seems missing??
env_cluster <- unique(df_predictors$cluster)
veg_cluster <- unique(df_stems$cluster)

sort(env_cluster)
sort(veg_cluster)

setdiff(env_cluster, veg_cluster)

pairs(sum_stems   ~ tmp + tmp_z + distance, df_fin)

# some values are missing predictors: eg. 18_178, 18_184

#5. test model ---------------------------------

# test glm with neg bin family
library(MASS)
global.negbin <- glm.nb(sum_stems ~ tmp  + tmp_z +  spei + prec + prcp_z + disturbance_severity + distance + manag + elevation, data = df_fin, na.action = "na.fail") # + tmp_z 

dredge(global.negbin)


global.negbin2 <- glm.nb(sum_stems ~ tmp   +  prec  + spei + disturbance_severity + distance, data = df_fin, na.action = "na.fail") # + tmp_z 

dredge(global.negbin2)

vif(global.negbin2)

plot(allEffects(global.negbin2))


m.nb <- glm.nb(sum_stems ~ tmp  +  distance + manag, data = df_fin, na.action = "na.fail") # + tmp_z 

m2 <- glm.nb(sum_stems ~ tmp   +  prec  + spei + disturbance_severity + distance + manag, data = df_fin, na.action = "na.fail") # + tmp_z 

m3 <- glm.nb(sum_stems ~ tmp_z  +  distance, data = df_fin, na.action = "na.fail") # + tmp_z 

simulationOutput <- simulateResiduals(fittedModel = m2, plot = T)


cor(df_fin$tmp, df_fin$tmp_z)

AIC(m1, m2)
# zero inflated model; hurdle model:
library(pscl)
# m.zip <- zeroinfl(sum_stems ~ tmp   +  prec  + disturbance_severity + distance, #+ tmp_z 
#                   dist = "poisson", data = df_fin)

m.zinb <- zeroinfl(sum_stems ~ tmp   +  spei + prec  + disturbance_severity + distance, #+ tmp_z
                   dist = "negbin", data = df_fin)

# m.hurdle1 is the best one for now!
m.hurdle1 <- hurdle(sum_stems ~ tmp   +  spei + prec  + disturbance_severity + distance, #+ tmp_z
                   dist = "negbin", data = df_fin)

m.hurdle2 <- hurdle(sum_stems ~ tmp   +  spei + prec  + disturbance_severity + distance + manag, #+ tmp_z
                   dist = "negbin", data = df_fin)

m.hurdle3 <- hurdle(sum_stems ~  spei + disturbance_severity + distance + manag, #+ tmp_z
                    dist = "negbin", data = df_fin)

m.hurdle4 <- hurdle(sum_stems ~  spei + disturbance_severity + distance, #+ tmp_z
                    dist = "negbin", data = df_fin)

# add back temperature and precipitation, as they are from 2021-2023, spei is from 2018-2020 - the best!
m.hurdle5 <- hurdle(sum_stems ~  tmp + prec + spei + disturbance_severity + distance + manag, #+ tmp_z
                    dist = "negbin", data = df_fin)


# remove TMP as correlated with spei
m.hurdle6 <- hurdle(sum_stems ~  prec + spei + disturbance_severity + distance + manag, #+ tmp_z
                    dist = "negbin", data = df_fin)

# add soil depth
m.hurdle7 <- hurdle(sum_stems ~  tmp + prec + spei + disturbance_severity + distance + manag + depth_extract, #+ tmp_z
                    dist = "negbin", data = df_fin)

# add soil depth, remove again tmp and prec as correlated with spei
m.hurdle8 <- hurdle(sum_stems ~  spei + disturbance_severity + distance + manag + depth_extract, #+ tmp_z
                    dist = "negbin", data = df_fin)



# for several groups:
kruskal.test(sum_stems ~ manag, data = df_fin)


# Perform the Dunn test
dunnTest <- dunn.test(df_fin$sum_stems, df_fin$manag, method="bonferroni")
dunnTest



AIC(m2, m.hurdle1, m.hurdle2, m.hurdle3, m.hurdle8, m.hurdle5, m.hurdle6, m.hurdle7)

summary(m.hurdle7)
r2(m2)

# check model assumptions
plot(resid(m.hurdle1) ~ fitted(m.hurdle1))

cor(df_fin$tmp, df_fin$spei)
cor(df_fin$prec, df_fin$spei)


# go fo negative binomial, zero inf

# filtered SPEI INf values

# visualize teh hurdle model: https://library.virginia.edu/data/articles/getting-started-with-hurdle-models
sum(predict(m.hurdle1, type = "prob")[,1])  # number of zeros in teh predictoed data

# First 5 expected counts
predict(m.hurdle1, type = "response")[1:5]


#library(countreg)

#rootogram(m.hurdle1) # fit up to count xx



# make plot for hurdle5

m.hurdle5

summary(m.hurdle5)

# Function to determine if a value is significant
is_significant <- function(significance) {
  ifelse(significance == " ", "Non-significant", "Significant")
}

# Prepare data for count model
count_model_summary <- data.frame(
  Predictor = c("tmp", "prec", "spei", "disturbance_severity", "distance", "managUnmanaged"),
  Estimate = c(0.1120775, 0.0013370, 0.3777615, 0.6380436, 0.0007405, 0.1164662),
  StdError = c(0.0377988, 0.0002272, 0.3064617, 0.1849542, 0.0002089, 0.0703305),
  Significance = c("**", "***", "", "***", "***", ".")
)
count_model_summary$SigColor <- sapply(count_model_summary$Significance, is_significant)

# Plot for count model
p_count <- ggplot(count_model_summary, aes(x = Predictor, y = Estimate, color = SigColor)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - StdError, ymax = Estimate + StdError), width = 0.1) +
  scale_color_manual(values = c("Non-significant" = "grey", "Significant" = "black")) +
  geom_hline(yintercept = 0, lty = 'dashed', col = 'grey') +
  theme_classic() +
  theme(legend.position = 'none') + 
  labs(title = "Count Model Coefficients",
       y = "Estimate",
       x = "",
       color = "Significance")

# Prepare data for zero hurdle model
zero_hurdle_model_summary <- data.frame(
  Predictor = c("tmp", "prec", "spei", "disturbance_severity", "distance", "managUnmanaged"),
  Estimate = c(0.429226, 0.002572, 5.470773, -0.397865, 0.001769, -0.146980),
  StdError = c(0.194297, 0.001051, 1.794975, 0.765920, 0.001345, 0.273274),
  Significance = c("*", "*", "**", "", "", "")
)
zero_hurdle_model_summary$SigColor <- sapply(zero_hurdle_model_summary$Significance, is_significant)

# Plot for zero hurdle model
p_zero_hurdle <- ggplot(zero_hurdle_model_summary, aes(x = Predictor, y = Estimate, color = SigColor)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - StdError, ymax = Estimate + StdError), width = 0.1) +
  scale_color_manual(values = c("Non-significant" = "grey", "Significant" = "black")) +
  theme_classic() +
  theme(legend.position = 'none') +
  geom_hline(yintercept = 0, lty = 'dashed', col = 'grey') +
  labs(title = "Zero Hurdle Model Coefficient",
       y = "Estimate",
       x = "",
       color = "Significance")

# Display the plots


ggarrange(p_count, p_zero_hurdle, nrow = 2)


# 6. management effect: compare manag vs unmanaged stem densities, richness---------------------

create_plot <- function(df, x_var, y_var, y_label, colors) {
  df <- df %>% 
    dplyr::filter(cluster %in% cluster_id_keep) 
  
  ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) + 
    stat_summary(fun = "median", geom = "bar") +
    stat_summary(fun.data = function(z) {
      list(y = median(z), ymin = quantile(z, 0.25), ymax = quantile(z, 0.75))
    }, geom = "errorbar", width = .2) +
    ylab(y_label) +
    xlab('') +
    scale_fill_manual(values = colors) +
    theme_classic() + 
    theme(legend.position = 'none')
}

# Define colors for 'manag' categories
# Assuming 'manag' has two distinct categories
# Replace with actual category names and desired colors
colors <- c("Managed" = "darkolivegreen", "Unmanaged" = "darkolivegreen1")

# Using the function
p.richness <- create_plot(df_richness, "manag", "richness", "Richness", colors)
p.density <- create_plot(df_stems, "manag", "sum_stems", "Stems [ha]", colors)
p.vert <- create_plot(df_vert, "manag", "n_layers", "# Vertical Layers", colors)
p.IVI.max <- create_plot(df_IVI_max, "manag", "rIVI", "Dominance (rIVI [%])", colors)

# Arrange plots
ggarrange(p.density, p.vert, p.IVI.max, p.richness, ncol = 4)


# test significance - for 2 groups only
wilcox.test(sum_stems ~ manag, data = df_stems)
# test significance - for 2 groups only
wilcox.test(n_layers ~ manag, data = df_vert)

# test significance - for 2 groups only
wilcox.test(rIVI ~ manag, data = df_IVI_max)

# test significance - for 2 groups only
wilcox.test(richness ~ manag, data = df_richness)




# Country effect ----------------------------------------------------------


df_richness <- df_richness %>% 
  dplyr::filter(cluster %in% cluster_id_keep) 


p.country.richness <- ggplot(data = df_richness,
       aes(x = country, y = richness,fill = country)) + 
  # Add bars as medians
  stat_summary(fun = "median", 
               geom = "bar") +
  stat_summary(
    data = df_richness,
    mapping = aes(x = country, y = richness),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) +
  ylab('Richness [#]') + 
  theme_classic() +
  theme(legend.position = 'none')
  #coord_cartesian(ylim=c(60, 67)) # ylim=c(59,66)



stem_dens_ha_cluster_sum <- stem_dens_ha_cluster_sum %>% 
  dplyr::filter(cluster %in% cluster_id_keep) 


df_stems_country <- stem_dens_ha_cluster_sum %>% 
  ungroup(.) %>% 
  group_by(cluster, manag, country) %>% # sum stems per luster, manag and country
  dplyr::reframe(sum_stems = sum(total_stems)) 

p.country.density <-  df_stems_country %>% 
  ggplot(aes(x = country, y = sum_stems,fill = country )) +
  # Add bars as medians
  stat_summary(fun = "median", 
               geom = "bar") +
  stat_summary(
    data = df_stems_country,
    mapping = aes(x = country, y = sum_stems),
    fun.min = function(z) { quantile(z,0.25) },
    fun.max = function(z) { quantile(z,0.75) },
    fun = median,
    geom  = 'errorbar',
    width = .2) +
  ylab('Stems [ha]') + 
  theme_classic() +
  theme(legend.position = 'none')

(p.country.density)
#coord_cartesian(ylim=c(60, 67)) # ylim=c(59,66)

# test significance - for 2 groups only
wilcox.test(sum_stems ~ manag, data = df_stems)


# for several groups:
kruskal.test(sum_stems ~ country, data = df_stems_country)


# Perform the Dunn test
dunnTest <- dunn.test(df_stems_country$sum_stems, df_stems_country$manag, method="bonferroni")
dunnTest
  

ggarrange(p.country.density, p.country.richness)



#### Stem density vs IVI -----------------------------------------------------

# Merge structure & composition into single space 

df_org_space <- 
  df_stems %>% 
  left_join(df_IVI_max, by = join_by(manag, cluster)) %>% 
  right_join(dplyr::select(df_richness, -country), by = join_by(manag, cluster)) %>% 
  inner_join(dplyr::select(df_vert, -sum_stems, -country), by = join_by(manag, cluster)) %>% 
  dplyr::filter(cluster %in% cluster_id_keep ) #%>% 
  #mutate(n_layers = case_when(is.na(n_layers)~0))

View(df_org_space)
# get summary for the scatter plot, one point is one country & management
df_summary <- df_org_space %>%
  group_by(manag) %>% # country, 
  summarize(
    rich_mean = median(rIVI, na.rm = TRUE),
    rich_sd = sd(rIVI, na.rm = TRUE),
    rich_min = quantile(rIVI, 0.25, na.rm = TRUE), # rich_mean -rich_sd, #
    rich_max = quantile(rIVI, 0.75, na.rm = TRUE), # rich_mean +rich_sd,
    dens_mean = median(sum_stems   , na.rm = TRUE),
    dens_sd   = sd(sum_stems  , na.rm = TRUE),
    dens_min = quantile(sum_stems  , 0.25, na.rm = TRUE), # dens_mean - dens_sd, #
    dens_max = quantile(sum_stems  , 0.75, na.rm = TRUE), # dens_mean + dens_sd, #
    .groups = 'drop'
  )

(df_summary)


df_org_space %>% 
  ggplot(aes(x = rIVI,
             y = sum_stems ,
             color = factor(manag)))+ 
  geom_point() +
  facet_wrap(.~country)



# prepare plot with raster

library(MASS)

# Function to calculate levels
calculate_levels <- function(density, percentages) {
  sorted_density <- sort(density, decreasing = TRUE)
  cum_density <- cumsum(sorted_density) / sum(sorted_density)
  sapply(percentages, function(p) {
    sorted_density[max(which(cum_density <= p))]
  })
}

# Calculate 2D density
d <- kde2d(df_org_space$rIVI, df_org_space$sum_stems, n = 500)
density_values <- d$z

# Set densities outside the range of your data to zero
density_values[density_values < 1e-6] <- 0

# Calculate the levels for specified percentages
levels <- calculate_levels(as.vector(density_values), c(0.50, 0.75, 0.90, 0.99, 1))

# Prepare data for ggplot
plot_data <- expand.grid(rIVI = d$x, sum_stems = d$y)
plot_data$density <- as.vector(density_values)

# Use cut to create factor levels, including one for zero density
plot_data$level <- cut(plot_data$density, breaks = c(-Inf, levels, Inf), labels = FALSE, include.lowest = TRUE)

# Define colors (from red to yellow, plus white for zero density)

library(RColorBrewer)
blue_colors <- brewer.pal(5, "Greens")

# Add 'white' at the beginning
my_colors <- c("white", blue_colors)





p.org.space.stems.IVI <- ggplot(plot_data) +
  geom_raster(aes(x = rIVI, y = sum_stems, fill = factor(level)), alpha = 0.8) +
  scale_fill_manual(values = my_colors,
                    labels = c("", rev(c("50%", "75%", "90%", "99%", "100%"))),
                    name = "Density") +
   ylab('Stem density') +
  xlab(bquote('Species dominance (rIVI)')) +
  theme_minimal() +
  #facet_wrap(.~country, scales = 'free') + 
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

 
# for vertical vs richness -------------------------------------------
# Calculate 2D density
d <- kde2d(df_org_space$richness, df_org_space$n_layers, n = 500)
density_values <- d$z

# Set densities outside the range of your data to zero
density_values[density_values < 1e-6] <- 0

# Calculate the levels for specified percentages
levels <- calculate_levels(as.vector(density_values), c(0.50, 0.75, 0.90, 0.99, 1))

# Prepare data for ggplot
plot_data <- expand.grid(richness = d$x, n_layers = d$y)
plot_data$density <- as.vector(density_values)

# Use cut to create factor levels, including one for zero density
plot_data$level <- cut(plot_data$density, breaks = c(-Inf, levels, Inf), labels = FALSE, include.lowest = TRUE)

# Define colors (from red to yellow, plus white for zero density)

library(RColorBrewer)
blue_colors <- brewer.pal(5, "Greens")

# Add 'white' at the beginning
my_colors <- c("white", blue_colors)





p.org.sp.rich.vert <- ggplot(plot_data) +
  geom_raster(aes(x = richness, y = n_layers, fill = factor(level)), alpha = 0.8) +
  scale_fill_manual(values = my_colors,
                    labels = c("", rev(c("50%", "75%", "90%", "99%", "100%"))),
                    name = "Density") +
  ylab('Vert. layers [#]') +
  xlab(bquote('Richness [#]')) +
  theme_minimal() +
  #facet_wrap(.~country, scales = 'free') + 
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

ggarrange(p.org.space.stems.IVI, p.org.sp.rich.vert,  ncol = 2, align = 'hv', common.legend = T, 
          legend = 'right' )

#           