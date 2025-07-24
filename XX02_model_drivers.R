# 5. Drivers ---------------------------------

# Model

gc()


library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(mgcv)
library(gratia)
library(ggeffects)
library(sjPlot)
library(corrplot)
library(spdep)


public_dir <- here("outData", "public")

# Load cleaned input data
df_fin <- fread(file.path(public_dir, "data", "plot_level_predictors_clean.csv"))
df_stem_species_class <- fread(file.path(public_dir, "data", "plot_level_stem_density_species_by_class.csv"))


# Variables
# total number of plots
n_total_plots = length(unique(df_fin$plot)) # 849

source('00_my_functions.R')


# select main variables as predictors 
predictor_vars_sub <- c(
  
  # over 2018-2023
  "spei1",
  "spei12",
  "tmp", 
  "prcp", 
  "tmp_z", 
  "prcp_z", 
  
  # during 2018-2020
  "drought_tmp",
  "drought_prcp",
  "drought_spei1",
  "drought_spei12",
  
  # disturbance chars
  "distance_edge", 
  "disturbance_severity", # from RS 
  
  # plot info
  "sand_extract",
  "clay_extract", 
  "depth_extract", 
  "av.nitro",
  
  # seasonality: CV - over year
  "cv_t2m",
  "cv_tp" #,
)


# Subset the data 

df_stem_regeneration <- df_fin %>% 
  dplyr::select(all_of(c("plot", 
                         "stem_regeneration",
                         predictor_vars_sub,
                         "country_pooled", 
                         "x", "y"))) %>% 
  mutate(country_pooled = factor(country_pooled))


summary(df_stem_regeneration)


## Pre-diagnostics --------------------------------
### Univariate models to find teh best predictors  -----------------------------

# Define response and predictor variables
response_var <- "stem_regeneration" # Replace with your actual response variable name
predictor_vars <- predictor_vars_sub

# Create an empty data frame to store results
AIC_results_univariate <- data.frame(
  Predictor = character(),
  AIC = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each predictor variable
for (predictor in predictor_vars) {
  # Create the formula
  formula <- as.formula(paste(response_var, "~ s(", predictor, ", k = 3)"))
  
  # Fit the GAM model with Tweedie distribution
  model <- gam(formula, family = tw(), data = df_stem_regeneration )
  
  # Extract AIC
  aic <- AIC(model)
  
  # Store the results
  AIC_results_univariate <- rbind(AIC_results_univariate, data.frame(Predictor = predictor, AIC = aic))
}

# Sort results by AIC (lower is better)
AIC_results_univariate <- AIC_results_univariate[order(AIC_results_univariate$AIC), ]

# Display the results
View(AIC_results_univariate)

sjPlot::tab_df(AIC_results_univariate,
               show.rownames = FALSE,
               file = file.path(public_dir, "tables", "table_univariate_GAM_AIC_comparison.doc"),
               digits = 1) 


### Spearman correlation coefficient ---------------------------------------------

# check correlation between predictors
library(corrplot)

# Select the relevant predictors from your data frame
predictors <- df_stem_regeneration %>%
  # dplyr::select(where(is.numeric))
  dplyr::select(all_of(predictor_vars_sub))

# Calculate the correlation matrix
correlation_matrix <- cor(predictors, use = "complete.obs", method = "spearman")

# Rename variables in the correlation matrix
correlation_matrix_renamed <- correlation_matrix

# Rename rows
rownames(correlation_matrix_renamed) <- rownames(correlation_matrix_renamed) %>%
  stringr::str_replace("^cv_t2m$", "cv_tmp") %>%
  stringr::str_replace("^cv_tp$", "cv_prcp")

# Rename columns
colnames(correlation_matrix_renamed) <- colnames(correlation_matrix_renamed) %>%
  stringr::str_replace("^cv_t2m$", "cv_tmp") %>%
  stringr::str_replace("^cv_tp$", "cv_prcp") %>% 
  stringr::str_replace("_extract$", "_content")


# Export with high resolution and clearer layout
png(file.path(public_dir, "figs", "fig_correlation_matrix_predictors.png"), 
    width = 2000, height = 2000, res = 300)

corrplot(correlation_matrix_renamed, 
         method = "color",
         type = "upper",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.col = "black",         # Label color
         tl.cex = 0.8,             # Smaller label font
         tl.srt = 45,              # Diagonal tilt (more readable)
         #addCoef.col = NA,         # Omit numbers for clarity
         # number.cex = 0.6,         # If you later add coefficients
         mar = c(0, 0, 1, 0)       # Reduce margins
)

dev.off()

## GAM(M)s------------

m_rnd <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(prcp,tmp, k = 5 ) +
    s(country_pooled, bs = "re") +
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration
)

m_rnd_te <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    #s(distance_edge, k = 5) +
    #s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    te(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(country_pooled, bs = "re") +
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration
)
m_rnd_ti <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(country_pooled, bs = "re") +
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration
)

m_rnd_ti_fixed <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5 ) +
    ti(prcp,tmp, k = 5 ) +
    s(x, y),                                 # Spatial autocorrelation
  family = tw(),
  method = 'REML',
  select = TRUE,
  data = df_stem_regeneration
)

AIC(m_rnd_ti_fixed, m_rnd_ti)

### Diagnostics  ----------------

mgcv::summary.gam(m_rnd_ti_fixed)
gratia::appraise(m_rnd_ti_fixed)


### Save model --------------------------------------------------
fin.m.reg.density <- m_rnd_ti_fixed

save(fin.m.reg.density, df_stem_regeneration,
     file = file.path(public_dir, "model", "model_stem_regeneration_density.RData"))



#### Relative contribution of factors  ------------------------------------
### drop one variable at time:  

# Fit full model
m_full <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5) +
    ti(prcp, tmp, k = 5) +
    s(x, y),
  family = tw(),
  method = "REML", select = TRUE,
  data = df_stem_regeneration
)
# Store logLik and AIC
loglik_full <- logLik(m_full)
aic_full <- AIC(m_full)

# Build and compare models one by one
m_drop_prcp <- gam(stem_regeneration ~ 
                     s(tmp, k = 5) + s(distance_edge, k = 5) +
                     s(disturbance_severity, k = 5) + s(clay_extract, k = 5) +
                     s(av.nitro, k = 5) + ti(disturbance_severity, distance_edge, k = 5) +
                     ti(prcp, tmp, k = 5) + s(x, y),
                   family = tw(), method = "REML",  select = TRUE, data = df_stem_regeneration)

m_drop_tmp <- gam(stem_regeneration ~ 
                    s(prcp, k = 5) + s(distance_edge, k = 5) +
                    s(disturbance_severity, k = 5) + s(clay_extract, k = 5) +
                    s(av.nitro, k = 5) + ti(disturbance_severity, distance_edge, k = 5) +
                    ti(prcp, tmp, k = 5) + s(x, y),
                  family = tw(), method = "REML",  select = TRUE,data = df_stem_regeneration)

m_drop_dist_edge <- gam(stem_regeneration ~ 
                          s(prcp, k = 5) + s(tmp, k = 5) +
                          s(disturbance_severity, k = 5) + s(clay_extract, k = 5) +
                          s(av.nitro, k = 5) + ti(disturbance_severity, distance_edge, k = 5) +
                          ti(prcp, tmp, k = 5) + s(x, y),
                        family = tw(), method = "REML",  select = TRUE, data = df_stem_regeneration)

m_drop_dist_sev <- gam(stem_regeneration ~ 
                         s(prcp, k = 5) + s(tmp, k = 5) +
                         s(distance_edge, k = 5) + s(clay_extract, k = 5) +
                         s(av.nitro, k = 5) + ti(prcp, tmp, k = 5) + s(x, y),
                       family = tw(), method = "REML",  select = TRUE, data = df_stem_regeneration)

m_drop_prcp_tmp <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5) +
    # ti(prcp, tmp, k = 5) +
    s(x, y),
  family = tw(),
  method = "REML",
  select = TRUE,
  data = df_stem_regeneration
)

# Compare models
results <- data.frame(
  model = c("m_rnd_ti_fixed", "drop_prcp", "drop_tmp", "drop_dist_edge", "drop_dist_sev", 'm_drop_prcp_tmp'),
  AIC = c(AIC(m_rnd_ti_fixed), AIC(m_drop_prcp), AIC(m_drop_tmp), AIC(m_drop_dist_edge), AIC(m_drop_dist_sev), AIC(m_drop_prcp_tmp)),
  logLik = c(logLik(m_rnd_ti_fixed), logLik(m_drop_prcp), logLik(m_drop_tmp), logLik(m_drop_dist_edge), logLik(m_drop_dist_sev), logLik(m_drop_prcp_tmp)),
  R2_adj = c(summary(m_rnd_ti_fixed)$r.sq, summary(m_drop_prcp)$r.sq, summary(m_drop_tmp)$r.sq,
             summary(m_drop_dist_edge)$r.sq, summary(m_drop_dist_sev)$r.sq, summary(m_drop_prcp_tmp)$r.sq),
  deviance_expl = c(summary(m_rnd_ti_fixed)$dev.expl, summary(m_drop_prcp)$dev.expl, summary(m_drop_tmp)$dev.expl,
                    summary(m_drop_dist_edge)$dev.expl, summary(m_drop_dist_sev)$dev.expl, summary(m_drop_prcp_tmp)$dev.expl)
)

results
# 
# > results
# model      AIC    logLik     R2_adj deviance_expl
# 1  m_rnd_ti_fixed 16069.08 -8013.875 0.09310789    0.10887421
# 2       drop_prcp 16072.64 -8011.886 0.09795147    0.11237287
# 3        drop_tmp 16070.49 -8013.018 0.09464378    0.11037273
# 4  drop_dist_edge 16073.63 -8015.650 0.09262081    0.10577277
# 5   drop_dist_sev 16076.56 -8019.264 0.09183355    0.09942951
# 6 m_drop_prcp_tmp 16070.75 -8015.014 0.09329267    0.10687028

### Low vs high severity plots ----------------------------

# Calculate the percentage of plots with >90% severity
prop_high_severity <- mean(df_stem_regeneration$disturbance_severity > 0.9, na.rm = TRUE) * 100
prop_high_severity

# Subset: Only plots with disturbance severity < 90%
df_low_severity <- df_stem_regeneration %>%
  dplyr::filter(disturbance_severity < 0.7)

nrow(df_low_severity)

# Refit the same GAM on this subset
m_lowsev <- gam(
  stem_regeneration ~ 
    s(prcp, k = 5) + s(tmp, k = 5) +
    s(distance_edge, k = 5) +
    s(disturbance_severity, k = 5) +
    s(clay_extract, k = 5) +
    s(av.nitro, k = 5) +
    ti(disturbance_severity, distance_edge, k = 5) +
    ti(prcp, tmp, k = 5) +
    s(x, y),
  family = tw(),
  method = "REML",
  select = TRUE,
  data = df_low_severity
)
plot(m_lowsev, pages = 1)
plot(m_rnd_ti_fixed, pages = 1)


summary(m_lowsev)
summary(m_rnd_ti_fixed)


#low vs full severity plots

# 1. Predict from the full model (across full severity range)
pred_full <- ggpredict(m_rnd_ti_fixed, terms = "disturbance_severity [0:1 by=0.01]") # 
pred_full$group <- "Full model (n = 849)"

# 2. Predict from the low-severity model (limited to 0–0.7 range)
pred_low <- ggpredict(m_lowsev, terms = "disturbance_severity [0:0.7 by=0.01]") #  
pred_low$group <- "Low-severity subset (n = 133)"

# 3. Combine predictions
pred_combined <- bind_rows(as.data.frame(pred_full), as.data.frame(pred_low)) %>%
  mutate(disturbance_severity_pct = x * 100)

# 4. Plot both prediction lines + confidence bands + raw data
p_low_high_severity <- 
  ggplot(pred_combined, aes(x = disturbance_severity_pct, y = predicted/1000, 
                            color = group,
                            linetype = group)) +
  geom_ribbon(aes(ymin = conf.low/1000, 
                  ymax = conf.high/1000, fill = group), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Full model (n = 849)" = "steelblue",
                                "Low-severity subset (n = 133)" = "darkorange"),
                     guide = guide_legend(nrow = 2)) +
  scale_fill_manual(values = c("Full model (n = 849)" = "steelblue",
                               "Low-severity subset (n = 133)" = "darkorange"),
                    guide = guide_legend(nrow = 2)) +
  scale_linetype_manual(values = c("Full model (n = 849)" = "solid",
                                   "Low-severity subset (n = 133)" = "dashed"),
                        guide = guide_legend(nrow = 2)) +
  
  labs(
    x = "Disturbance severity [%]",
    y = expression("Stem regeneration density [*1000 n ha"^-1*"]"),
    color = "", fill = "", linetype = "",
  ) +
  theme_classic2(base_size = 8) +
  theme(legend.position = "bottom")


p_low_high_severity

# Save the last plot
ggsave(p_low_high_severity,
       filename = file.path(public_dir, "figs", "figS_stem_density_by_disturbance_severity_subset_vs_full.png"),
       width = 4,
       height = 4,
       units = "in",
       dpi = 300  # High resolution for publication
)






### Export final drivers model: only fixed effectd  -----------------------------------------------------------------

# Identify random effects using the model's "smooth" component
smooth_terms <- summary(fin.m.reg.density)$s.table
# Identify random effects using the model's "smooth" component
smooth_terms_lowsev <- summary(m_lowsev)$s.table


# Extract the smooth terms labels and check which ones are random effects
random_effects_labels <- rownames(smooth_terms)[str_detect(rownames(smooth_terms), "country_pooled|clim_grid")]
random_effects_labels_lowsev <- rownames(smooth_terms_lowsev)[str_detect(rownames(smooth_terms_lowsev), "country_pooled|clim_grid")]


# Create a function to automatically label the terms in the summary output
create_labels <- function(term) {
  if (term %in% random_effects_labels) {
    return(paste(term, "(Random Effect)"))
  } else {
    return(term)
  }
}

# Apply the labeling function
pred_labels <- sapply(rownames(smooth_terms), create_labels)
pred_labels_lowsev <- sapply(rownames(smooth_terms_lowsev), create_labels)

# Save tab_model for full model
sjPlot::tab_model(
  fin.m.reg.density,
  show.re.var = TRUE,
  pred.labels = c("Intercept", pred_labels),
  dv.labels = paste0("Explained Deviance: ", round(100 * summary(fin.m.reg.density)$dev.expl, 2), "%"),
  file = file.path(public_dir, "tables", "table_model_stem_regeneration_density_full.doc")
)

# Save tab_model for low-severity subset model
sjPlot::tab_model(
  m_lowsev,
  show.re.var = TRUE,
  pred.labels = c("Intercept", pred_labels_lowsev),
  dv.labels = paste0("Explained Deviance: ", round(100 * summary(m_lowsev)$dev.expl, 2), "%"),
  file = file.path(public_dir, "tables", "table_model_stem_regeneration_density_low_severity.doc")
)

#### with random effects 

# Identify random effects using the model's "smooth" component
smooth_terms <- summary(m_rnd_ti)$s.table

# Extract the smooth terms labels and check which ones are random effects
random_effects_labels <- rownames(smooth_terms)[str_detect(rownames(smooth_terms), "country_pooled|clim_grid")]

# Create a function to automatically label the terms in the summary output
create_labels <- function(term) {
  if (term %in% random_effects_labels) {
    return(paste(term, "(Random Effect)"))
  } else {
    return(term)
  }
}

# Apply the labeling function
pred_labels <- sapply(rownames(smooth_terms), create_labels)

# Display the tab_model with automatic labeling
sjPlot::tab_model(
  m_rnd_ti,
  show.re.var = TRUE,
  pred.labels = c("Intercept", pred_labels),
  dv.labels = paste0("Explained Deviance: ", round(100 * summary(m_rnd_ti)$dev.expl, 2), "%"),
  file = file.path(public_dir, "tables", "table_model_stem_regeneration_density_full_with_random_effects.doc")
)




### Plot: Drivers  ---------------------------------------------------------------------------

y_lab = expression("Stem density [1000 n ha"^{-1}*"]")
y_lab_reg = expression("Reg. stem density [1000 n ha"^{-1}*"]")

m <- fin.m.reg.density #m_int_sev_edge_full_te_comb # m_int_res_edge_full_te_comb #  m_int_res_edge_full_te

# Generate predictions using ggpredict
summary(m)

# Interaction 1: Precipitation and Temperature
pred_prcp       <- ggpredict(m, terms = c("prcp"))
pred_tmp        <- ggpredict(m, terms = c("tmp"))
pred_tmp_prcp   <- ggpredict(m, terms = c("prcp", "tmp [8,10]"))
pred_dist_edge  <- ggpredict(m, terms = c("distance_edge[50:250]"))
pred_dist_sever <- ggpredict(m, terms = c("disturbance_severity"))
pred_clay       <- ggpredict(m, terms = c("clay_extract"))

# Convert all ggpredict objects to data.frames
pred_prcp        <- as.data.frame(pred_prcp)
pred_tmp         <- as.data.frame(pred_tmp)
pred_tmp_prcp    <- as.data.frame(pred_tmp_prcp)
pred_dist_edge   <- as.data.frame(pred_dist_edge)
pred_dist_sever  <- as.data.frame(pred_dist_sever)
pred_clay        <- as.data.frame(pred_clay)


my_theme_drivers <- theme(
  axis.title = element_text(size = 8),
  plot.title = element_text(hjust = 0.5, size = 8),       # Title size
  #axis.title.y = element_blank(),                   # Axis title size
  axis.text = element_text(size = 8),                    # Axis text size
  legend.key.size = unit(0.5, "cm"),                     # Legend key size
  legend.text = element_text(size = 8),                  # Legend text size
  legend.title = element_text(size = 8),                 # Legend title size
  legend.position = c(0.05, 0.9),
  legend.justification = c(0.05, 0.9)
)

#my_colors_interaction <- c("grey90", "red") 
my_colors_interaction <- c("#FDAE61", "#A50026") 
my_color_main_effects <- "grey90" # "#006837"    

# Plot precipitation
p.prcp <- ggplot(pred_prcp, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), alpha = 0.4,fill = my_color_main_effects) +
  labs(x = "Precipitation [mm]", 
       y = y_lab_reg, 
       title = "p<0.001") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p.prcp

# Plot precipitation
p.tmp <- ggplot(pred_tmp, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), alpha = 0.4,fill = my_color_main_effects) +
  labs(x = "Temperature [°C]",  
       y = y_lab_reg, 
       title = "p=0.004") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p.tmp



# Plot the first interaction
p1 <- 
  ggplot(pred_tmp_prcp, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, aes(color = group) ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000, fill = group), 
              alpha = 0.2, color = NA) +
  scale_color_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  scale_fill_manual(values = my_colors_interaction, name = "Temperature [°C]") +
  # theme_classic() +
  labs(x = "Precipitation [mm]", 
       y =y_lab_reg,  
       title = "p=0.006", 
       # linetype =  "Temperature [°C]"
  ) +
  my_theme_drivers

p1

# Plot distance to edge
p2 <- ggplot(pred_dist_edge, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), 
              alpha = 0.3,fill = my_color_main_effects) +
  labs(x = "Distance to edge [m]",  
       y = y_lab_reg, 
       title = "p=0.023") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p2


# 
p3 <- ggplot(pred_dist_sever, aes(x = x*100, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000), fill =my_color_main_effects, alpha = 0.3, color = NA) +
  theme_classic() +
  # ylim(0,15) +
  labs(x = "Disturbance severity [%]", 
       y = "", title = "p<0.001") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p3


# PlotCaly
p4 <- ggplot(pred_clay, aes(x = x, y = predicted/1000)) +
  geom_line(linewidth = 1, color = my_color_main_effects ) +
  geom_ribbon(aes(ymin = conf.low/1000, ymax = conf.high/1000),  
              fill = my_color_main_effects, alpha = 0.3, color = NA) +
  theme_classic() +
  scale_y_continuous(breaks = seq(5, 15, 5)) +
  #ylim(0,20) +
  labs(x = "Clay content [%]", y = "", title = "p=0.003") +
  my_theme_drivers + 
  theme(legend.position = 'none')
p4




p_combined_int_no_points <- ggarrange(p1, p4, p2, p3,
                                      labels = c("[a]","[b]", "[c]","[d]"), 
                                      align = 'hv',
                                      font.label = list(size = 8, face = "plain")) # Specify plain font style)

p_combined_int_no_points

# Save the combined plot
ggsave(
  filename = file.path(public_dir, "figs", "Fig2.png"),
  plot = p_combined_int_no_points,
  width = 5.5,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)


#### all predictors for supplement  ------------------------

# update y lables plotting
p.tmp <- p.tmp + labs(y = "")
p1    <- p1 + labs(y = "")
p2    <- p2 + labs(y = "")
p3    <- p3 + labs(y = "")
p4    <- p4 + labs(y = y_lab_reg )


p_combined_int_no_points_supplem <- ggarrange(p.prcp, p.tmp, p1, p4, p2, p3,
                                              labels = c("[a]","[b]", "[c]","[d]", "[e]","[f]"), 
                                              ncol = 3, nrow = 2,
                                              align = 'hv',
                                              font.label = list(size = 8, face = "plain")) # Specify plain font style)

p_combined_int_no_points_supplem

# Save the combined plot
ggsave(
  filename = file.path(public_dir, "figs", "figS_regeneration_by_drivers_supplement.png"),
  plot = p_combined_int_no_points_supplem,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)




