
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


# cluster analysis: ------------------------------------------------------------
# for cliamtic and environmental condistions to be on the same page with Kilian's clusters

# env-climate condition) -------------------
# Subset the relevant columns
data_subset_clim <- df_indicators[, c("tmp", "prcp", "tmp_z", "prcp_z", "spei", "sand_extract", "clay_extract", "depth_extract", "av.nitro")]
#data_subset_str  <- df_indicators[, c("dominant_species", "rIVI", "richness", "management_intensity", "stem_density", "n_vertical", "distance_edge", "disturbance_severity")]

#
# Standardize the data
data_scaled_clim <- scale(data_subset_clim)
#data_scaled_str  <- scale(data_subset_str)
#data_scaled_str  <- df_indicators$dominant_species

# Find which number of clusters is teh best
# Perform K-means clustering for different values of k
set.seed(3)
max_clusters <- 10
sil_width <- numeric(max_clusters)
for (k in 2:max_clusters) {
  kmeans_result_clim <- kmeans(data_scaled_clim, centers = k, nstart = 25)
  sil <- silhouette(kmeans_result_clim$cluster, dist(data_scaled_clim))
  sil_width[k] <- mean(sil[, 3])
}

# Plot Silhouette width for different values of k
plot(1:max_clusters, sil_width, type = "b", xlab = "Number of clusters", ylab = "Average Silhouette width", main = "Silhouette Analysis for K-means Clustering")

# Determine the optimal number of clusters
optimal_k_clim <- 3  # from Kilian's study

# Perform K-means clustering with the optimal number of clusters
set.seed(3)
kmeans_result <- kmeans(data_scaled_clim, centers = optimal_k_clim, nstart = 25)

# Add cluster assignments to the original data
df_indicators$clim_cluster <- kmeans_result$cluster

# Perform PCA for visualization - use Pc1 and Pc1
pca_result_clim <- prcomp(data_scaled_clim)

# plot PCA results with the most important variables: 
biplot(pca_result_clim, main = "PCA Biplot")

pairs(data_scaled_clim)

library(GGally) # for pairwise comparions of the variables 1:1, also calculates the correlation and sinificance

# Convert to data frame if it's not already
data_scaled_clim_df <- as.data.frame(data_scaled_clim)

# Create a ggpairs plot
ggpairs(data_scaled_clim_df)

#### Plotting ENV_CLIM cluster analysis --------------------------------------------------
# Plot the PCA results with clusters
# plot(pca_result$x[, 1:2], col = kmeans_result$cluster, pch = 20, main = "K-means Clustering Result (PCA)", xlab = "PC1", ylab = "PC2")
# points(kmeans_result$centers %*% pca_result$rotation[, 1:2], col = 1:optimal_k, pch = 8, cex = 2)
# 
# # Add legend
# legend("topright", legend = paste("Cluster", 1:optimal_k), col = 1:optimal_k, pch = 8)
# 
# 
# # color the PCa by the cluster numebrs:
# # Plot the PCA results with clusters
# plot(pca_result$x[, 1], pca_result$x[, 2], col = kmeans_result$cluster, pch = 20, 
#      main = "K-means Clustering Result (PCA)", xlab = "PC1", ylab = "PC2")
# points(kmeans_result$centers %*% pca_result$rotation[, 1:2], col = 1:optimal_k, pch = 8, cex = 2)
# 
# 
# # Visualize the clustering result
# plot(data_scaled, col = kmeans_result$cluster, pch = 20, main = "K-means Clustering Result")
# points(kmeans_result$centers, col = 1:optimal_k, pch = 8, cex = 2)
# legend("topright", legend = paste("Cluster", 1:optimal_k), col = 1:optimal_k, pch = 8)
# 
# 
# 
# 
# # Visualize the clustering result
# plot(x = data_scaled[, "tmp_z"], y = data_scaled[, "spei"], col = kmeans_result$cluster, pch = 20, main = "K-means Clustering Result")
# legend("topright", legend = paste("Cluster", 1:optimal_k), col = 1:optimal_k, pch = 8)
# 
# # Extract cluster centers for the specified x and y axes
# centers_scaled <- scale(kmeans_result$centers, center = attr(data_scaled, "scaled:center"), scale = attr(data_scaled, "scaled:scale"))
# 
# # Plot the cluster centers on the same axes
# points(centers_scaled[, "tmp_z"], centers_scaled[, "spei"], col = 1:optimal_k, pch = 8, cex = 2)
# # Add legend
# legend("topright", legend = paste("Cluster", 1:optimal_k), col = 1:optimal_k, pch = 8)
# 

# Silhouette plot for the chosen number of clusters
# silhouette_result <- silhouette(kmeans_result$cluster, dist(data_scaled))
# plot(silhouette_result, main = "Silhouette Plot", col = as.numeric(silhouette_result[, 1]))

# Cluster analysis for structural characteristics -------------------------------
# first split data in 3 atasets, run the cluste analysis for each one of them
# from Kilian: 4, 5, 2 clusters for 1.2.3 clim clusters;
# need to used k-prototypes
# Add cluster assignments to the original data
df_indicators$clim_cluster <- kmeans_result$cluster



# structural cluster analysis 
# Subset the data
# "dominant_species", 
data_subset_str  <- df_indicators[, c("dominant_species", "rIVI", "richness", "management_intensity", 
                                      "stem_density", "n_vertical", "distance_edge", "disturbance_severity",
                                      "clim_cluster")]

# Convert dominant_species to factor
data_subset_str$dominant_species <- as.factor(data_subset_str$dominant_species)

# Normalize the quantitative variables for clustering
data_normalized <- data_subset_str %>%
  mutate(across(c(rIVI, richness, management_intensity, stem_density, n_vertical, distance_edge, disturbance_severity), scale)) %>% 
  as.data.frame()

# Dummy code the categorical variable
data_dummied <- model.matrix(~dominant_species - 1, data = data_normalized)
data_dummied <- as.data.frame(data_dummied)
# Combine the dummy coded categorical variable with the normalized quantitative variables
data_for_clustering <- cbind(data_dummied, data_normalized[, -1])

# avoid dummy variables, hard to read
# split table in 3 clusters to run separately structurral cluster analysis:
df1 <- data_normalized[, -1] %>%  dplyr::filter(clim_cluster == 1)
df2 <- data_normalized[, -1] %>%  dplyr::filter(clim_cluster == 2)
df3 <- data_normalized[, -1] %>%  dplyr::filter(clim_cluster == 3)

# remove last columns so they are not part of teh PCA aalysis
df1 <- df1 %>% dplyr::select(-clim_cluster)
df2 <- df2 %>% dplyr::select(-clim_cluster)
df3 <- df3 %>% dplyr::select(-clim_cluster)




# Perform clustering using k-means (you can choose other methods like PAM if desired)
set.seed(123) # Setting seed for reproducibility
kmeans_result1 <- kmeans(df1, centers = 4, nstart = 25) # Adjust centers as needed
kmeans_result2 <- kmeans(df2, centers = 5, nstart = 25) # Adjust centers as needed
kmeans_result3 <- kmeans(df3, centers = 3, nstart = 25) # Adjust centers as needed

# Perform PCA for visualization - use Pc1 and Pc1
pca_result_str1 <- prcomp(df1)#prcomp(data_for_clustering)
pca_result_str2 <- prcomp(df2)#prcomp(data_for_clustering)
pca_result_str3 <- prcomp(df3)#prcomp(data_for_clustering)

# plot PCA results with the most important variables: 
#biplot(pca_result_str1, main = "PCA Biplot")


library(ggbiplot)


p1 <- ggbiplot(pca_result_str1,varname.color = "red", circle = TRUE) + ggtitle('1 = cold')
p2 <- ggbiplot(pca_result_str2,varname.color = "red", circle = TRUE) + ggtitle('2 = hot')
p3 <- ggbiplot(pca_result_str3,varname.color = "red", circle = TRUE) + ggtitle('3 = med')
#ggscreeplot(pca_result_str3)

ggarrange(p1, p2, p3, nrow = 3, ncol = 1)
p1
p2
p3



# Add cluster assignments to the original data
df1$str_cluster <- kmeans_result1$cluster
df2$str_cluster <- kmeans_result2$cluster
df3$str_cluster <- kmeans_result3$cluster

df <- rbind(df1, df2, df3)





# Get country indications ----------------------------------------------------------

# add country indication 
country_regions <- tribble(
  ~country, ~regions, ~country_abbr,
  "germany", "11, 12, 14, 18, 19, 20, 25", "DE",
  "poland", "17", "PL",
  "czech", "15, 26", "CZ",
  "austria", "13", "AT",
  "slovakia", "16", "SK",
  "slovenia", "23", "SI",
  "italy", "21", "IT",
  "switzerland", "22", "CH",
  "france", "24, 27", "FR"
) %>%
  separate_rows(regions, sep = ", ") %>%
  mutate(region = as.integer(regions))


# Extract the region number from the site column in df_indicators
df_indicators <- df_indicators %>%
  mutate(region = as.integer(substr(site, 1, 2)))

# Merge the country information with df_indicators
df_indicators_with_country <- df_indicators %>%
  left_join(country_regions, by = "region") %>%
  mutate(country_abbr = case_when(
    site %in% c("24_136", "24_137") ~ "BE",  # Belgium
    site %in% c("24_133", "24_134", "24_135") ~ "LU",  # Luxembourg
    TRUE ~ country_abbr
  )) %>% 
  mutate(country_pooled = case_when( country_abbr == "BE" ~ "FR",  # create pooled data for the eco analysis
                                    country_abbr == "LX" ~ "FR",
                                    TRUE~country_abbr)) %>% 
  mutate(clim_class = case_when(clim_cluster  == '1' ~ 'wet-cold',
                                clim_cluster  == '2'  ~ 'hot-dry',
                                clim_cluster  == '3'  ~ 'medium',
                                TRUE ~ NA_character_
  ))  
  

# investigate the clim clusters by country
df_indicators_with_country %>% 
  ggplot(aes(x = country_abbr,
             fill = factor(clim_class))) +
  geom_bar(position = 'dodge')

head(df_indicators)


# split the ENV-CLIM clusters into each group, as each one can have a different number of clusters
# data into 3 groups





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


# Make a scatter plot: how does the spei correlate with temperature?

df_indicators_with_country %>% 
  ggplot(aes(x = tmp, y = spei, color = clim_class )) + 
  geom_point()





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
