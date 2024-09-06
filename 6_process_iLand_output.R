
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

# Process input data -----------------------------------------------------------
# get simulated data
df_sim <- fread('outTable/df_simulated.csv')
#View(df_sim)

# get field data to compare with simulated ones
df_field     <- fread('outData/veg_density_DBH.csv')
df_indicators <- fread('outData/indicators_for_cluster_analysis.csv')

df_indicators <- df_indicators %>% 
  rename(prcp = prec) %>% 
  rename(site = cluster)


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



# filter field data for selected clusters (simulated landscapes) --------------
df_field_sub <- df_field %>% 
  rename(site = cluster) %>% 
  right_join(df_sites_clusters) %>% 
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)


# get plot for all of locations:
df_indicators <- df_indicators %>% 
 # rename(site = cluster) %>% 
  left_join(df_sites_clusters) %>%
 # mutate()
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)


df_indicators_sub <- df_indicators %>% 
#  rename(site = cluster) %>% 
  right_join(df_sites_clusters) %>% 
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)


df_sim <- df_sim %>%
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)




# Get summary stats for 3 clusters ----------------------------------------------------
# Calculate summary statistics grouped by 'clim_cluster_test'
summary_stats_clim_cluster_country <- 
  df_indicators_with_country %>%
  group_by(clim_class, country_abbr) %>%
  summarise(
    tmp_min = min(tmp, na.rm = TRUE),
    tmp_max = max(tmp, na.rm = TRUE),
    tmp_median = median(tmp, na.rm = TRUE),
    tmp_mean = mean(tmp, na.rm = TRUE),
    tmp_sd = sd(tmp, na.rm = TRUE),
    #tmp_iqr = IQR(tmp, na.rm = TRUE),
    prcp_min = min(prcp, na.rm = TRUE),
    prcp_max = max(prcp, na.rm = TRUE),
    prcp_median = median(prcp, na.rm = TRUE),
    prcp_mean = mean(prcp, na.rm = TRUE),
    prcp_sd = sd(prcp, na.rm = TRUE),
    spei_min = min(spei, na.rm = TRUE),
    spei_max = max(spei, na.rm = TRUE),
    spei_median = median(spei, na.rm = TRUE),
    spei_mean = mean(spei, na.rm = TRUE),
    spei_sd = sd(spei, na.rm = TRUE),
    #spei_iqr = IQR(spei, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -c(clim_class,country_abbr),
    names_to = c("variable", "stat"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>% 
  mutate(
      mean_sd = paste0(round(mean, 1), " ± ", round(sd, 1)),
      min = round(min, 1),
      max = round(max, 1),
      median = round(median, 1)
    ) %>%
    dplyr::select(clim_class, country_abbr, variable, mean_sd, median, min, max) # Reorder columns
  

  
sjPlot::tab_df(summary_stats_clim_cluster_country,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = F,
               file="outTable/clim_cluster_summary_country.doc",
               digits = 2) 


# Sumary stat per climate cluster
summary_stats_clim_cluster <- 
  df_indicators_with_country %>%
  group_by(clim_class) %>%
  summarise(
    tmp_min = min(tmp, na.rm = TRUE),
    tmp_max = max(tmp, na.rm = TRUE),
    tmp_median = median(tmp, na.rm = TRUE),
    tmp_mean = mean(tmp, na.rm = TRUE),
    tmp_sd = sd(tmp, na.rm = TRUE),
    #tmp_iqr = IQR(tmp, na.rm = TRUE),
    prcp_min = min(prcp, na.rm = TRUE),
    prcp_max = max(prcp, na.rm = TRUE),
    prcp_median = median(prcp, na.rm = TRUE),
    prcp_mean = mean(prcp, na.rm = TRUE),
    prcp_sd = sd(prcp, na.rm = TRUE),
    spei_min = min(spei, na.rm = TRUE),
    spei_max = max(spei, na.rm = TRUE),
    spei_median = median(spei, na.rm = TRUE),
    spei_mean = mean(spei, na.rm = TRUE),
    spei_sd = sd(spei, na.rm = TRUE),
    #spei_iqr = IQR(spei, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = -c(clim_class),
    names_to = c("variable", "stat"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>% 
  mutate(
    mean_sd = paste0(round(mean, 1), " ± ", round(sd, 1)),
    min = round(min, 1),
    max = round(max, 1),
    median = round(median, 1)
  ) %>%
  dplyr::select(clim_class,variable, mean_sd, median, min, max) # Reorder columns



sjPlot::tab_df(summary_stats_clim_cluster,
               #col.header = c(as.character(qntils), 'mean'),
               show.rownames = F,
               file="outTable/clim_cluster_summary.doc",
               digits = 1) 


# Make a scatter plot: how does the spei correlate with temperature -----------------------------------

# Calculate the mean and standard deviation for each cluster
mean_sd_df <- df_indicators_with_country %>%
  group_by(clim_class) %>%
  summarise(
    mean_tmp = mean(tmp, na.rm = TRUE),
    sd_tmp = sd(tmp, na.rm = TRUE),
    mean_spei = mean(spei, na.rm = TRUE),
    sd_spei = sd(spei, na.rm = TRUE)
  )



# Create the scatter plot with ellipses, mean points, and error bars
fig_spei_tmp_clusters <- ggplot(data = df_indicators_with_country, aes(x = tmp, y = spei, color = clim_class, fill = clim_class)) + 
  geom_point(size = 0.5) +
  stat_ellipse(aes(group = clim_class), type = "norm", alpha = 0.3, geom = "polygon") +
  geom_point(data = mean_sd_df, aes(x = mean_tmp, y = mean_spei), size = 2.8, shape = 16, color = 'white') +  # Add halo efferct
  
  geom_point(data = mean_sd_df, aes(x = mean_tmp, y = mean_spei), size = 1.5, shape = 16) +  # Mean points
  geom_errorbar(data = mean_sd_df, aes(x = mean_tmp, ymin = mean_spei - sd_spei, ymax = mean_spei + sd_spei, y = mean_spei), 
                width = 0.2) +  # Vertical error bars
  geom_errorbarh(data = mean_sd_df, aes(y = mean_spei, xmin = mean_tmp - sd_tmp, xmax = mean_tmp + sd_tmp, x = mean_tmp), 
                 height = 0.01) +  # Horizontal error bars
  theme_classic2() +
  labs(title = "",
       x = expression(Temperature ~ (degree*C)),
       y = "SPEI") +
  scale_fill_manual(values = c("red", "orange", "blue")) +  # Customize fill colors
  scale_color_manual(values = c("red", "orange", "blue"))   # Customize point colors




(fig_spei_tmp_clusters )
ggsave(filename = 'outFigs/fig_clim_clusters.png', plot = fig_spei_tmp_clusters, 
       width = 4, height = 3, dpi = 300, 
       bg = 'white')



### Make climate space for each country ---------------------------------------------

mean_sd_df_country <- df_indicators_with_country %>%
  group_by(clim_class, country_pooled) %>%
  summarise(
    mean_tmp = mean(tmp, na.rm = TRUE),
    sd_tmp = sd(tmp, na.rm = TRUE),
    mean_spei = mean(spei, na.rm = TRUE),
    sd_spei = sd(spei, na.rm = TRUE)
  )




fig_spei_tmp_clusters_country <- 
  ggplot(data = df_indicators_with_country, aes(x = tmp, y = spei, color = clim_class, fill = clim_class)) + 
  geom_point(size = 0.5) +
  stat_ellipse(aes(group = clim_class), type = "norm", alpha = 0.3, geom = "polygon") +
  geom_point(data = mean_sd_df_country, aes(x = mean_tmp, y = mean_spei), size = 2.8, shape = 16, color = 'white') +  # Add halo effect
  geom_point(data = mean_sd_df_country, aes(x = mean_tmp, y = mean_spei), size = 1.5, shape = 16) +  # Mean points
  geom_errorbar(data = mean_sd_df_country, aes(x = mean_tmp, ymin = mean_spei - sd_spei, ymax = mean_spei + sd_spei, y = mean_spei), 
                width = 0) +  # Vertical error bars
  geom_errorbarh(data = mean_sd_df_country, aes(y = mean_spei, xmin = mean_tmp - sd_tmp, xmax = mean_tmp + sd_tmp, x = mean_tmp), 
                 height = 0) +  # Horizontal error bars
  theme_classic2() +
  labs(title = "",
       x = expression(Temperature ~ (degree*C)),
       y = "SPEI") +
  scale_fill_manual(values = c("red", "orange", "blue")) +  # Customize fill colors
  scale_color_manual(values = c("red", "orange", "blue")) +  # Customize point colors
  facet_wrap(~ country_pooled)  # Add facet_wrap to create separate panels for each country

# Display and save the plot
print(fig_spei_tmp_clusters_country)
ggsave(filename = 'outFigs/fig_spei_tmp_clusters_country.png', plot = fig_spei_tmp_clusters_country, 
       width = 7, height = 7, dpi = 300, 
       bg = 'white')






# # overall stats
summary_stats_all <- df_indicators %>%
  ungroup(.) %>% 
 # group_by(clim_cluster_test) %>%
  summarise(
    tmp_min = min(tmp, na.rm = TRUE),
    tmp_max = max(tmp, na.rm = TRUE),
    tmp_median = median(tmp, na.rm = TRUE),
    tmp_mean = mean(tmp, na.rm = TRUE),
    tmp_sd = sd(tmp, na.rm = TRUE),
    tmp_iqr = IQR(tmp, na.rm = TRUE),
    prcp_min = min(prcp, na.rm = TRUE),
    prcp_max = max(prcp, na.rm = TRUE),
    prcp_median = median(prcp, na.rm = TRUE),
    prcp_mean = mean(prcp, na.rm = TRUE),
    prcp_sd = sd(prcp, na.rm = TRUE),
    prcp_iqr = IQR(prcp, na.rm = TRUE),
    spei_min = min(spei, na.rm = TRUE),
    spei_max = max(spei, na.rm = TRUE),
    spei_median = median(spei, na.rm = TRUE),
    spei_mean = mean(spei, na.rm = TRUE),
    spei_sd = sd(spei, na.rm = TRUE),
    spei_iqr = IQR(spei, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", "stat"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )

summary_stats_all

### Inspect clim drivers: For all sites  -------------------
# to understand which cluster is what:
# Create an empty list to store the plots
plot_list_all <- list()

# Loop through combinations and create plots
for (combination in plot_combinations) {
  p <- ggplot(df_indicators, aes_string(x = combination$x, y = combination$y)) +
    geom_point(color = 'grey', alpha = 0.3) + 
   # geom_point(data = df_indicators_sub, aes_string(x = combination$x, y = combination$y, color = clim_cluster), size = 4) +
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
 

# make a plot manually: tmp spei
p1 <- ggplot(df_indicators, aes(x = tmp, y = spei)) +
  geom_point(color = 'grey', alpha = 0.3) + 
  geom_text(data = df_indicators_sub, aes(x = tmp, y = spei, label = cluster, color = clim_cluster), 
            vjust = -0.1, size = 6) +
  geom_point(data = df_indicators_sub, aes(x = tmp, y = spei, color = clim_cluster), size = 3) +
  # labs(#title = paste(combination$x, "vs", combination$y),
  #   x = "tmp",
  #   y = "spei") +
  theme_bw() +
  theme(aspect.ratio = 1)  # Ensure the plot is square


p2 <- ggplot(df_indicators, aes(x = tmp_z, y = spei)) +
  geom_point(color = 'grey', alpha = 0.3) + 
  geom_text(data = df_indicators_sub, aes(x = tmp_z, y = spei, label = cluster, color = clim_cluster), 
            vjust = -0.1, size = 6) +
  geom_point(data = df_indicators_sub, aes(x = tmp_z, y = spei, color = clim_cluster), size = 3) +
   theme_bw() +
  theme(aspect.ratio = 1)  # Ensure the plot is square


p3 <- ggplot(df_indicators, aes(x = tmp_z, y = prcp_z)) +
  geom_point(color = 'grey', alpha = 0.3) + 
  geom_text(data = df_indicators_sub, aes(x = tmp_z, y = prcp_z, label = cluster, color = clim_cluster), 
            vjust = -0.1, size = 6) +
  geom_point(data = df_indicators_sub, aes(x = tmp_z, y = prcp_z, color = clim_cluster), size = 3) +
  theme_bw() +
  theme(aspect.ratio = 1)  # Ensure the plot is square

ggarrange(p1,p2,p3, common.legend = T, ncol = 3, nrow = 1)

# plot only subset: climatic clusters (12 points)
df_indicators_sub %>% 
  ggplot(aes(x = spei,
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

# Process simulated data --------------------------------------------------------
### Inspect simulated data: ------------------------------------------------------
# 12 clusters, 
# 3 climatic models
# 4 climate scenarios
# 2 seed added/no
unique(df_sim$cluster)  # 12, 3 climatic, 4 strcutural in each climate cluster
#[1] "1_1" "1_2" "1_3" "1_4" "2_1" "2_2" "2_3" "2_4" "2_5" "3_1" "3_2" "3_3"

unique(df_sim$species) # 24
#[1] "acps" "cabe" "potr" "saca" "quro" "acca" "frex" "ulgl" "lade" "tico" "piab" "fasy" "psme"
#[14] "soau" "pisy" "bepe" "abal" "alin" "algl" "soar" "coav" "acpl" "rops" "casa"


unique(df_sim$clim_model)
# [1] "ICHEC" "MPI"   "NCC"  

unique(df_sim$clim_scenario)
# [1] "HISTO" "RCP26" "RCP45" "RCP85"

# external seed: external seed present or not?? TRUE/FALSE
# years: 1-30
unique(df_sim$run_nr)  # 5 repetitions


# df_cl %>% 
#   dplyr::filter(category == "unclassified") %>% 
#   summary()




### Indicators from simulated ata : ------------------------------------------------------ 
# richness
# rIVI
# stem density
# vertical classes
# summarize across simulation repetition and across the clim model

##### Species richness ------------------------------------
df_richness <- df_sim %>% 
  group_by(year, cluster,  clim_scenario, ext_seed) %>%  #clim_modelclim_cluster, 
  summarize(richness = n_distinct(species))

df_richness %>% 
 # dplyr::filter(run_nr == 1 & clim_model == 'NCC')  %>% 
  ggplot() + 
  geom_line(aes(x = year,
                y = richness,
                color = cluster,
                group = cluster)) + 
  facet_wrap(clim_scenario~ext_seed)


##### Species importance value (relative, from rel_density and rel_BA) ------------------------------------
# get species importance value: from relative density, relative BA
# first calculate the total values per ha, then add it to original table to calculate teh rIVI based on relative dominance
df_sum_cluster <- df_sim %>% 
  group_by(year,  cluster, clim_scenario, ext_seed) %>%# clim_cluster clim_model ,,
  summarize(sum_stems = sum(count_ha, na.rm = T),
            sum_BA = sum(basal_area_m2, na.rm = T)) %>% 
  ungroup()


df_IVI <- df_sim %>% 
  # group by speies across the levels
  group_by(year, species, cluster, clim_scenario, ext_seed) %>% #clim_cluster,clim_model , 
  summarize(sp_dens = sum(count_ha),
            sp_BA   = sum(basal_area_m2)) %>% 
  ungroup(.) %>% 
    left_join(df_sum_cluster) %>% # merge sums per whole cluster
  mutate(rel_dens  = sp_dens/sum_stems,
         rel_BA  = sp_BA/sum_BA,
         rIVI = ( rel_dens +rel_BA)/2) %>% # relative species importance value
  # filter teh dominant species - highest rIVI
  group_by(year, cluster,  clim_scenario, ext_seed) %>% #clim_model ,
  dplyr::filter(rIVI == max(rIVI)) %>%
  sample_n(1) %>%  # Select a random row if there are ties
  rename(dominant_species      = species)
  
  


##### Structure: Vertical classes & stem density -------------------------------------------------------------
# inspeact vertical classes
df_structure <- df_sim %>% 
  group_by(year, cluster, clim_scenario, ext_seed) %>% #clim_model , 
  summarize(n_vertical = n_distinct(category),
            stem_density = sum(count_ha)) 


# change in number of vertical layers
ggplot(df_structure) + 
  geom_line(aes(x = year,
                y = n_vertical,
                color = cluster,
                group = cluster)) + 
  facet_wrap(clim_scenario~ext_seed)



# change in stem density
ggplot(df_structure) + 
  geom_line(aes(x = year,
                y = stem_density,
                color = cluster,
                group = cluster)) + 
  facet_wrap(clim_scenario~ext_seed)


# merge df indicators 
df_sim_indicators <- df_richness %>% 
  left_join(df_IVI, by = join_by(year, cluster, clim_scenario, ext_seed)) %>%
  left_join(df_structure, by = join_by(year, cluster, clim_scenario, ext_seed)) %>% 
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1))  #%>% # add indication of the strutural cluster (1,2,3,4)
 
  

  
##### evaluate initial state with my sites (only 12 sites! )   -------------------------
# filterr initial state: year == 0
df_sim_indicators0 <- df_sim_indicators %>% 
  dplyr::filter(year == 0 & clim_scenario == "HISTO")


###### merge field with simulated data in year 0 
df_compare <- df_indicators_sub %>% 
  left_join(df_sim_indicators0, by = "cluster", suffix = c("_field", "_simul")) %>% 
  dplyr::select(cluster, site, ext_seed, ends_with("_field"), ends_with("_simul")) %>% 
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1))  %>% # add indication of the strutural cluster (1,2,3,4)
 na.omit() # remove empty one
  
  

print(df_compare)


# Create a list of plots
plot_list <- list()

# Extract column names with "_field" and "_simul" suffixes
field_columns <- grep("_field$", names(df_compare), value = TRUE)
simul_columns <- gsub("_field$", "_simul", field_columns)

# Generate scatter plots for each pair of columns
for (i in seq_along(field_columns)) {
  p <- ggplot(df_compare, aes_string(x = field_columns[i], y = simul_columns[i])) +
    geom_point(aes(color = clim_cluster), size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +  # Add diagonal line
    labs(x = field_columns[i], y = simul_columns[i], 
         title =""  #paste(field_columns[i], "vs", simul_columns[i])
         ) +
    theme_bw()
  
  plot_list[[i]] <- p
}

# Print all plots
p1 <- plot_list[[1]]
p2 <- plot_list[[2]]
p3 <- plot_list[[3]]
p4 <- plot_list[[4]]
p5 <- plot_list[[5]]

windows()
ggarrange(p1, p2, p3, p4,p5, common.legend = TRUE , ncol = 3, nrow = 2)


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
df_indicators_sub <- 
  df_indicators %>% 
  dplyr::select(site, clim_cluster_test, dominant_species,rIVI, richness, stem_density, n_vertical) %>% 
  rename(clim_cluster = clim_cluster_test) %>% 
  mutate(year = 0,
         clim_scenario  = "HISTO") %>% 
  dplyr::select(year, site, clim_cluster, clim_scenario, dominant_species, rIVI, richness, stem_density, n_vertical) %>% 
  mutate(source = 'observed')

 
# merge tables into one
df_indi_merged <- rbind(df_indicators_sub,df_sim_indicators_avg_sub)

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
