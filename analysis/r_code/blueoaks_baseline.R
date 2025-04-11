source(here::here("analysis", "r_code", "helpers", "analysis_helpers.R"))
source(here::here("analysis", "r_code", "helpers", "plot_helpers.R"))

nc <- parallelly::availableCores()
options(mc.cores = nc)

set.seed(1)

calib_start_lobos <- 1950
# recon_start <- 1800

w_djf <- c(1/3, 1/3, rep(0, 9), 1/3)

monthly_precip_lobos <- read_csv(file = here('data', 'prism_casestudy.csv'))

# do transformations to climate first
lobos_precip_df <- monthly_precip_lobos %>%
  filter(year < 2024) %>%
  mutate(log_precip = log(prcp + 0.001)) %>%
  dplyr::select(
    year,
    month,
    log_precip,
    precip = 'prcp'
  )


# aggregate to get our reconstruction target
clim_df_djf_lobos <- process_climate(lobos_precip_df,
                                      clim_name = 'precip',
                                      w_agg = w_djf, anom = F)


# take log at yearly level
clim_df_djf_lobos <- clim_df_djf_lobos %>%
  mutate(
    log_clim = log(clim + 0.001)
  )

# pull out mean and standard deviation
# so we can standardize only using the training period
logclim_training <- clim_df_djf_lobos %>%
  filter(year >= calib_start_lobos) %>%
  pull(log_clim)

mean_logtarget = mean(logclim_training)
sd_logtarget = sd(logclim_training)


# standardize
clim_df_djf_lobos_std <- clim_df_djf_lobos %>%
  mutate(
    clim = (log_clim - mean_logtarget) / sd_logtarget
  ) %>%
  select(-log_clim)

trw_dfs_lobos <- format_trw_lobos(start_year_calib = calib_start_lobos)

trw_dfs_mod_lobos <- process_trw(trw_dfs_lobos)


df_list_lobos_baseline <- stan_setup(clim_df_djf_lobos_std, trw_dfs_mod_lobos,
                            start_year_calib = calib_start_lobos, month_model = F)

df_list_lobos_baseline$mean_logdjf = mean_logtarget
df_list_lobos_baseline$sd_logdjf = sd_logtarget


# model fit, including a model discrepancy term
fit_lobos_baseline <- stan(file = here("analysis", "stan_code", "baseline_model_eta.stan"),
                           iter= 5000,
                           data = df_list_lobos_baseline,
                           init_r = 1,
                           seed = 123)

# saveRDS(fit_lobos_baseline, here('temp_output', 'blueoaks_baseline_eta.Rds'))


##################

# model evaluation

####################

recon_start_lobos_baseline <- min(df_list_lobos_baseline$start_years)
validation_start_lobos_baseline <- max(recon_start_lobos_baseline, min(lobos_precip_df$year))

# extract reconstruction summary.
clim_raw_lobos_baseline <- fit_lobos_baseline %>%
  spread_draws(clim_raw[i]) %>%
  summarise_draws()

# get instrumental record for validation
djf_df_lobos <- process_climate(lobos_precip_df,
                              clim_name = 'precip',
                              w_agg = w_djf)

# dataframe with reconstruction and instrumental record (one row per year)
baseline_raw_lobos_df <- clim_raw_lobos_baseline %>%
  ungroup() %>%
  mutate(year = recon_start_lobos_baseline + i - 1,
         clim_mean_pred = mean
  ) %>%
  left_join(djf_df_lobos, by = 'year')

# plot reconstruction
fig_recon(baseline_raw_lobos_df,
          validation_start_lobos_baseline,
          calib_start_lobos,
          fig_only = T) +
  labs(title = str_glue('blue oaks'),
       y = 'DJF mean precip (mm/month)')



## repeat above, but in standardized log space --------

clim_mis_baseline <- fit_lobos_baseline %>%
  spread_draws(clim_yearly[i]) %>%
  summarise_draws()

baseline_lobos_df <- clim_mis_baseline %>%
  ungroup() %>%
  mutate(year = recon_start_lobos_baseline + i - 1,
         clim_mean_pred = mean
  ) %>%
  left_join(clim_df_djf_lobos_std, by = 'year') %>%
  rename(truth = clim)


fig_recon(baseline_lobos_df, validation_start_lobos_baseline, calib_start_lobos, truth_name = 'truth',
          display_years = c(1800, 1955), fig_only = T) +
  labs(title = 'blue oaks, baseline')

