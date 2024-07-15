
# Process: 
# Get indicators from simulated data
# - read iLand simulated data - simulated by Kilian on 29/06/2024
# - get vertical classes
# - calculate indicators form simulated data
# - average values between no seed and seed scenario - complete later by sensitivity analysis
# - evaluate they cchange change over time, how thtable they are given the clim cluster

# Process my data:
# - run k means clustering to classify clim_cluster
# - compare my indicators with 

# pre-analysis
# - inspect the development of 12 iLand landscapes


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
#library(factoextra)
#library(NbClust)


# get simulated data
df_sim <- fread('outTable/df_simulated.csv')
#View(df_sim)

# get field data to compare with simulated ones
df_field     <- fread('outData/veg_density_DBH.csv')
df_indicators <- fread('outData/indicators_for_cluster_analysis.csv')

df_indicators <- df_indicators %>% 
  rename(prcp = prec)
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

# Create the data frame with the given pairs
df_sites_clusters <- data.frame(
  site = c("23_132", "26_134", "15_133", "17_104", "22_101", "12_151",
           "24_146", "20_116", "12_117", "11_145", "19_160", "25_150"),
  cluster = c("1_1", "1_2", "1_3", "1_4", "2_1", "2_2",
              "2_3", "2_4", "2_5", "3_1", "3_2", "3_3")
)


# Cluster analysis based on environmental condistion(climate) -------------------
# Subset the relevant columns
data_subset <- df_indicators[, c("tmp", "prcp", "tmp_z", "prcp_z", "spei", "sand_extract", "clay_extract", "depth_extract", "av.nitro")]

# Standardize the data
data_scaled <- scale(data_subset)

# Perform K-means clustering for different values of k
set.seed(3)
max_clusters <- 10
sil_width <- numeric(max_clusters)
for (k in 2:max_clusters) {
  kmeans_result <- kmeans(data_scaled, centers = k, nstart = 25)
  sil <- silhouette(kmeans_result$cluster, dist(data_scaled))
  sil_width[k] <- mean(sil[, 3])
}

# Plot Silhouette width for different values of k
plot(1:max_clusters, sil_width, type = "b", xlab = "Number of clusters", ylab = "Average Silhouette width", main = "Silhouette Analysis for K-means Clustering")

# Determine the optimal number of clusters
optimal_k <- 3  # from Kilian's study


# Perform K-means clustering with the optimal number of clusters
set.seed(3)
kmeans_result <- kmeans(data_scaled, centers = optimal_k, nstart = 25)

# Add cluster assignments to the original data
df_indicators$clim_cluster_test <- kmeans_result$cluster

# Perform PCA for visualization - use Pc1 and Pc1
pca_result <- prcomp(data_scaled)

# plot PCA results with the most important variables: 
biplot(pca_result, main = "PCA Biplot")


# Plot the PCA results with clusters
plot(pca_result$x[, 1:2], col = kmeans_result$cluster, pch = 20, main = "K-means Clustering Result (PCA)", xlab = "PC1", ylab = "PC2")
points(kmeans_result$centers %*% pca_result$rotation[, 1:2], col = 1:optimal_k, pch = 8, cex = 2)

# Add legend
legend("topright", legend = paste("Cluster", 1:optimal_k), col = 1:optimal_k, pch = 8)


# color the PCa by the cluster numebrs:
# Plot the PCA results with clusters
plot(pca_result$x[, 1], pca_result$x[, 2], col = kmeans_result$cluster, pch = 20, 
     main = "K-means Clustering Result (PCA)", xlab = "PC1", ylab = "PC2")
points(kmeans_result$centers %*% pca_result$rotation[, 1:2], col = 1:optimal_k, pch = 8, cex = 2)


# Visualize the clustering result
plot(data_scaled, col = kmeans_result$cluster, pch = 20, main = "K-means Clustering Result")
points(kmeans_result$centers, col = 1:optimal_k, pch = 8, cex = 2)
legend("topright", legend = paste("Cluster", 1:optimal_k), col = 1:optimal_k, pch = 8)




# Visualize the clustering result
plot(x = data_scaled[, "tmp_z"], y = data_scaled[, "spei"], col = kmeans_result$cluster, pch = 20, main = "K-means Clustering Result")
legend("topright", legend = paste("Cluster", 1:optimal_k), col = 1:optimal_k, pch = 8)

# Extract cluster centers for the specified x and y axes
centers_scaled <- scale(kmeans_result$centers, center = attr(data_scaled, "scaled:center"), scale = attr(data_scaled, "scaled:scale"))

# Plot the cluster centers on the same axes
points(centers_scaled[, "tmp_z"], centers_scaled[, "spei"], col = 1:optimal_k, pch = 8, cex = 2)
# Add legend
legend("topright", legend = paste("Cluster", 1:optimal_k), col = 1:optimal_k, pch = 8)


# Silhouette plot for the chosen number of clusters
silhouette_result <- silhouette(kmeans_result$cluster, dist(data_scaled))
plot(silhouette_result, main = "Silhouette Plot", col = as.numeric(silhouette_result[, 1]))




# filter field data for selected clusters (simulated landscapes) --------------
df_field_sub <- df_field %>% 
  rename(site = cluster) %>% 
  right_join(df_sites_clusters) %>% 
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)


# get plot for all of locations:
df_indicators <- df_indicators %>% 
  rename(site = cluster) %>% 
  left_join(df_sites_clusters) %>%
 # mutate()
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)


df_indicators_sub <- df_indicators %>% 
  rename(site = cluster) %>% 
  right_join(df_sites_clusters) %>% 
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)


df_sim <- df_sim %>%
  mutate(clim_cluster = str_sub(cluster, 1, 1),  # add indication of the climatic cluster (1,2,3)
         str_cluster = str_sub(cluster, -1, -1))  # add indication of the strutural cluster (1,2,3,4,5)




# Plots: only for cluster data: analyse clim clusters tmp and prcp? ----------------------
#summary(df_indicators_sub) 


# List of combinations to plot
plot_combinations <- list(
  list(x = "tmp", y = "prcp"),
  list(x = "tmp_z", y = "prcp_z"),
  list(x = "tmp", y = "spei"),
  list(x = "tmp_z", y = "spei"),
  list(x = "prcp", y = "spei"),
  list(x = "prcp_z", y = "spei")
)

# Create an empty list to store the plots
plot_list_sub <- list()

# Loop through combinations and create plots
for (combination in plot_combinations) {
  p <- df_indicators_sub %>%
    ggplot(aes_string(x = combination$x, y = combination$y)) +
    geom_point(aes(color = clim_cluster)) + 
    geom_smooth(method = "loess", se = FALSE) +
    labs(title = paste(combination$x, "vs", combination$y),
         x = combination$x,
         y = combination$y) +
    theme_bw() +
    theme(aspect.ratio = 1)  # Ensure the plot is square
  
  plot_list_sub[[paste(combination$x, combination$y, sep = "_vs_")]] <- p
}

# Arrange and print all plots using ggarrange
ggarrange(plotlist = plot_list_sub, ncol = 2, nrow = 3, common.legend = T)


### For all sites  -------------------

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


# Check simulation data: ------------------------------------------------------
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




# calculate my variables on site level: ------------------------------------------------------ 
# richness
# rIVI
# stem density
# vertical classes
# summarize across simulation repetition and across the clim model

### Species richness ------------------------------------
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


### Species importance value (relative, from rel_density and rel_BA) ------------------------------------
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
  
  


# Structure: Vertical classes & stem density -------------------------------------------------------------
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
 
  

  
  
# filterr initial state: year == 0
df_sim_indicators0 <- df_sim_indicators %>% 
  dplyr::filter(year == 0 & clim_scenario == "HISTO")


# merge field with simulated data in year 0 ----------------------------------

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
