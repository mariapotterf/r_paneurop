
# Evaluate simulated forest growth development -----


# Process: 
# Get indicators from simulated data
# - read iLand simulated data - simulated by Kilian on 06/11/2024

# Calculate indicators from simulated data
# temporal develpment
#- evaluate they cchange change over time, how thtable they are given the clim cluster

# simulated data:
# clim_model: 3
# clim_scenario: 4 
# seeds range: 8
# repetition: 5
# landscapes: 12


rm(list = ls())             # Remove all objects from global environment
graphics.off()              # Close all open graphics devices
gc()                        # Run garbage collection to free up memory




library(data.table)  # fread()
library(dplyr)       # mutate(), filter(), group_by(), summarize()
library(ggplot2)     # ggplot(), geom_*
library(tidyr)       # Not used in this snippet, can be removed
library(stringr)     # str_extract(), str_remove()
library(ggpubr)      # For theme_classic2()

# Input data -----------------------------------------------------------------------------
public_dir <- here("outData", "public")

df_sim_indicators <- fread(file.path(public_dir, "data", "df_sim_indicators.csv"))

source(file.path(public_dir, "code", "00_my_functions.R"))

# Vars 

reg_colors_short <- c(
  "Del." = "#A50026",   # reddish
  "Int."   = "#FDAE61",   # yellowish
  "Adv."= "#006837"    # green
)


# prepare summary table
df_sens_plot <- df_sim_indicators %>% 
  # filter final years
  dplyr::filter(year %in% 25:30) %>%
  mutate(
    seed_comb = ifelse(is.na(seed_comb), "No seed", as.character(seed_comb)),
    seed_level_num = as.integer(str_extract(seed_comb, "-?\\d+")),
    seed_comb = factor(seed_comb, levels = c("No seed", paste0("seed_", 
                                                               sort(unique(seed_level_num))))),
    adv_delayed = recode_factor(adv_delayed,
                                "Delayed" = "Del.",
                                "Intermediate" = "Int.",
                                "Advanced" = "Adv."),
    adv_delayed = factor(adv_delayed, levels = c("Del.", "Int.", "Adv."))
  )

# calculate medians and IQR
group_summary <- df_sens_plot %>%
  dplyr::group_by(adv_delayed, seed_comb) %>%
  dplyr::summarize(
    q25 = quantile(stem_density / 1000, 0.25, na.rm = TRUE),
    median_density = median(stem_density / 1000, na.rm = TRUE),
    q75 = quantile(stem_density / 1000, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# make a plot
p_simulated_stem_sensitivity <- group_summary %>%
  ggplot(aes(y = seed_comb, x = median_density,
             xmin = q25, xmax = q75,
             color = adv_delayed)) +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.3) +
  geom_vline(data = group_summary %>%
               group_by(adv_delayed) %>%
               summarize(median_density = median(median_density), .groups = "drop"),
             aes(xintercept = median_density, color = adv_delayed),
             linetype = "dashed", show.legend = FALSE) +
  scale_color_manual(values = reg_colors_short) +
  theme_classic2() +
  scale_x_continuous(limits = c(0, 13)) +
  scale_y_discrete(labels = ~ str_remove(., "seed_")) +
  labs(y = expression("Change in seeds availability [%]"),
       x = expression("Stem density [" * 1000 *~ n~ha^{-1} * "]"),
       color = "") +
  theme(
    legend.position = 'bottom',
    panel.border = element_rect(color = "black", linewidth = 0.7, fill = NA),
    text = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.title = element_text(size = 8)
  )

p_simulated_stem_sensitivity


ggsave(
  filename = file.path(public_dir, "figs", "figS_simulated_sensitivity.png"), 
  plot = p_simulated_stem_sensitivity, 
  width = 3, 
  height = 4, 
  dpi = 300, 
  bg = 'white'
)
