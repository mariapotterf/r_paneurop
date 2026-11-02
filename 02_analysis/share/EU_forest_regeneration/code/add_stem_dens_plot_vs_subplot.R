### REGENERATION ONLY (saplings + juveniles): stems per subplot and plot, overall and per country ----
# Input: sub_plots_counts_DBH.csv  (one row per subplot x VegType x species, zeros included)
#   ID = subplot, cluster = plot, country, VegType (Saplings / Juveniles / Mature), count
# Resolution: subplot = 4 m2  -> 1 stem = 2500 n/ha
#             plot    = n_subplots x 4 m2 (normally 5 x 4 = 20 m2) -> 1 stem = 500 n/ha

library(data.table)
library(dplyr)
library(tidyr)

counts_file <- "outData/sub_plots_counts_DBH.csv"     # adjust path

subplot_area_m2  <- 4
fixed_5_subplots <- FALSE   # FALSE: plot area = subplots actually in the file x 4 m2
# TRUE : plot area = 20 m2 for every plot
veg_types <- c("Saplings", "Juveniles")               # regeneration only (mature trees excluded)

raw <- fread(counts_file) %>%
  filter(VegType %in% veg_types)

# 1) SUBPLOT level: stems per subplot (empty subplots = 0 rows already in the data) ----
subplot_df <- raw %>%
  group_by(country, plot = cluster, subplot = ID) %>%
  summarise(stems = sum(count, na.rm = TRUE), .groups = "drop") %>%
  mutate(stems_ha = stems / subplot_area_m2 * 10000)          # 4 m2 -> ha

# 2) PLOT level: sum of its subplots, scaled by the area of the plot ------------------
plot_df <- subplot_df %>%
  group_by(country, plot) %>%
  summarise(stems = sum(stems), n_subplots = n(), .groups = "drop") %>%
  mutate(area_m2  = if (fixed_5_subplots) 5 * subplot_area_m2 else n_subplots * subplot_area_m2,
         stems_ha = stems / area_m2 * 10000)                  # 20 m2 (or 16 m2) -> ha

# checks
n_distinct(plot_df$plot)        # expect 849 plots
nrow(subplot_df)                # expect 4221 subplots
sum(subplot_df$stems)           # total regeneration stems (12085)
table(plot_df$n_subplots)       # 24 plots have only 4 subplots
table(plot_df$country)          # countries in the file

# 3) summaries: overall + per country, plot + subplot, counts + n/ha ------------------
long <- bind_rows(
  plot_df    %>% transmute(country, level = "plot",    count = stems, per_ha = stems_ha),
  subplot_df %>% transmute(country, level = "subplot", count = stems, per_ha = stems_ha)
)

long_all <- bind_rows(long, long %>% mutate(country = "All countries"))

stem_summary <- long_all %>%
  pivot_longer(c(count, per_ha), names_to = "unit", values_to = "value") %>%
  group_by(country, level, unit) %>%
  summarise(
    n              = n(),
    mean           = mean(value),
    sd             = sd(value),
    median         = median(value),
    q25            = quantile(value, 0.25, names = FALSE),
    q75            = quantile(value, 0.75, names = FALSE),
    min            = min(value),
    max            = max(value),
    share_zero     = mean(value == 0) * 100,          # % empty plots / subplots
    mean_nonzero   = mean(value[value > 0]),
    median_nonzero = median(value[value > 0]),
    .groups = "drop"
  ) %>%
  mutate(country = factor(country,
                          levels = c("All countries", sort(unique(long$country))))) %>%
  arrange(country, level, unit)

print(stem_summary, n = Inf, width = Inf)

# compact view: mean and median only
stem_summary %>%
  select(country, level, unit, n, mean, median) %>%
  print(n = Inf)

# Table: means and medians, counts and n/ha, plot and subplot level ------------------
table_regen <- stem_summary %>%
  select(country, level, unit, mean, median) %>%
  pivot_wider(names_from = c(level, unit), values_from = c(mean, median),
              names_glue = "{level}_{unit}_{.value}") %>%
  select(country,
         plot_count_mean,    plot_count_median,    plot_per_ha_mean,    plot_per_ha_median,
         subplot_count_mean, subplot_count_median, subplot_per_ha_mean, subplot_per_ha_median) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))

print(table_regen, width = Inf)

# write.csv(table_regen, "outTable/stem_table_regeneration.csv", row.names = FALSE)
# write.csv(stem_summary, "outTable/stem_summary_plot_subplot_country.csv", row.names = FALSE)


# Boxplots: distribution per country (jitter dots drawn over the boxes) ---------------
# long_all has all countries + "All countries", plot and subplot level, counts and n/ha
library(ggpubr)
# Boxplots: distribution per country; plot and subplot level as groups ---------------
# jitter dots over the boxes; red diamond + number = MEDIAN of each plot / subplot group

ylim_ha    <- c(0, 40000)   # y-axis limits (zoom only, no data removed)
ylim_count <- c(0, 50)

plot_long <- long_all %>%
  pivot_longer(c(count, per_ha), names_to = "unit", values_to = "value") %>%
  mutate(
    country = factor(country, levels = c("All countries", sort(unique(long$country)))),
    level   = factor(level, levels = c("plot", "subplot"),
                     labels = c("Plot\n(20 m²)", "Subplot\n(4 m²)"))
  )

# median per country x level x unit (computed explicitly, so each group has its own median)
mean_df <- plot_long %>%
  group_by(unit, country, level) %>%
  summarise(mean = mean(value), .groups = "drop") %>%
  mutate(lab = if_else(unit == "per_ha",
                       format(round(mean), big.mark = ",", trim = TRUE),
                       as.character(round(mean, 1))))

# how many points fall above the y limit (hidden by the zoom)?
plot_long %>%
  filter(country == "All countries") %>%
  mutate(lim = if_else(unit == "per_ha", ylim_ha[2], ylim_count[2])) %>%
  group_by(unit, level) %>%
  summarise(pct_above_limit = mean(value > lim) * 100, .groups = "drop")

level_cols <- c("Plot\n(20 m²)" = "#74c364", "Subplot\n(4 m²)" = "#f4a736")
dodge_w    <- 0.75

make_box <- function(unit_name, ylab, ylim) {
  ggplot(filter(plot_long, unit == unit_name),
         aes(x = country, y = value, fill = level)) +
    geom_boxplot(position = position_dodge(width = dodge_w), width = 0.65,
                 outlier.shape = NA, colour = "grey20", alpha = 0.5, linewidth = 0.3) +
    # geom_jitter(position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0,
    #                                             dodge.width = dodge_w),
    #             size = 0.3, alpha = 0.25, colour = "grey20") +           # dots over the box
    geom_point(data = filter(mean_df, unit == unit_name),
               aes(y = mean, group = level),
               position = position_dodge(width = dodge_w),
               shape = 23, size = 1.8, #fill = "red", 
               colour = "black", 
               stroke = 0.3) +
    # geom_text(data = filter(median_df, unit == unit_name),
    #           aes(y = median, group = level, label = lab),
    #           position = position_dodge(width = dodge_w),
    #           vjust = -0.9, size = 2, colour = "black") +
    # scale_fill_manual(values = level_cols, name = NULL) +
    coord_cartesian(ylim = ylim) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
}

p_box_ha    <- make_box("per_ha", "Stem density [n/ha]",                  ylim_ha)
p_box_count <- make_box("count",  "Stems [count per plot / subplot]",     ylim_count)

p_box_countries <- ggarrange(p_box_ha, p_box_count, 
                             ncol =1,
                             nrow = 2,
                             common.legend = TRUE, legend = "bottom",
                             labels = c("[a]", "[b]"),
                             font.label = list(size = 10, face = "plain"))
p_box_countries

# ggsave(file.path(public_dir, "figs", "Fig_S_stems_boxplot_countries.png"),
#        p_box_countries, width = 7.5, height = 3.5, dpi = 300, units = "in")
