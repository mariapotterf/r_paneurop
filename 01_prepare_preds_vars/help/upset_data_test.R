#### upset data test ---------------------------------------------------------

dd <- data.frame(site = c(1,2,2,3,3,3,4,5,5,6,7,7,7),
                 vert = c('m', 
                          's', 'j',
                          'j','s','m',
                          's',
                          'j','s',
                          's',
                          'm','j','s'))
# Step 1: Create a binary presence/absence matrix for each site
dd_wide <- dd %>%
  pivot_wider(names_from = vert, values_from = vert, 
              values_fn = length, values_fill = 0) %>%
  mutate(m = ifelse(m > 0, 1, 0),
         j = ifelse(j > 0, 1, 0),
         s = ifelse(s > 0, 1, 0))

# Step 2: Select only the columns with presence/absence data
upset_data <- dd_wide %>% dplyr::select(m, j, s) %>% 
  as.data.frame()

# Step 3: Create the UpSet plot
upset(upset_data, sets = c("m", "j", "s"), order.by = "freq")





