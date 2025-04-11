source(here::here("analysis", "r_code", "helpers", "analysis_helpers.R"))
source(here::here("analysis", "r_code", "helpers", "plot_helpers.R"))

nc <- parallelly::availableCores()
options(mc.cores = nc)


set.seed(1)

calib_start_ppe <- 1950
w_djf <- c(1/3, 1/3, rep(0, 9), 1/3)

# monthly climate model output (first run CESM1 last mil)
sim_df <- read_csv(here('data', 'sim_cesm1_mil_data.csv')) %>%
  filter(year >= start_year_reconstruct) %>%
  filter(year != 2004) # 2004 only has data up to sept, so remove for simplicity to have only years of full data

# do transformations to climate first
sim_precip_df <- sim_df %>%
  mutate(log_precip = log(precip_mmmonth + 0.001),
         month = str_to_title(name)) %>%
  select(
    year,
    month,
    log_precip,
    precip = precip_mmmonth
  )

clim_df_ppe <- process_climate(sim_precip_df,
                               clim_name = 'log_precip',
                               anom = T)


# ordered jan - december
month_clim_ppe <- sim_precip_df %>%
  mutate(month = as_factor(month)) %>%
  group_by(month) %>%
  summarise(
    clim_mean = mean(log_precip),
    clim_sd = sd(log_precip)
  )

month_means_ppe <- month_clim_ppe %>%
  pull(clim_mean)
month_sds_ppe <- month_clim_ppe %>%
  pull(clim_sd)



snr_ppe <- 0.75
trw_dfs_ppe <- format_trw_linear(clim_df_ppe,
                                      n_trees = 10,
                                      snr = snr_ppe,
                                      start_year_calib = calib_start_ppe,
                                      w_true = w_djf, beta0 = 2, beta1 = 3, beta_age = -1)

trw_mod_dfs_ppe <- process_trw(trw_dfs_ppe)

halfyear_prior <- c(10,10,rep(1, 6), rep(10, 4))
unif_prior <- rep(1,12)
df_list_ppe <- stan_setup(clim_df_ppe, trw_mod_dfs_ppe,
                               start_year_calib = calib_start_ppe,
                          w_param = unif_prior, month_means = month_means_ppe,
                          month_sds = month_sds_ppe)


# init_r helps avoid non-posdef matrix initialization
# this will take a while to run!
fit_ppe_linear <- stan(file = here("analysis", "stan_code", "proposed_model.stan"),
                 iter= 5000,
                 data = df_list_ppe,
                init_r = 1,
                seed = 123#,
                #control = list(adapt_delta = 0.99)
                )

## code if you want to save this output!
# saveRDS(fit_ppe_linear, here('temp_output', 'fit_ppe_linear.Rds'))
# fit_ppe_linear <- readRDS(here('paper_output', 'fit_ppe_linear1950.Rds'))


fig_weight(fit_ppe_linear) +
  labs(title = str_glue('linear ppe with log(precip + 0.001), snr = {snr_ppe}'))


# reconstruction starts with first year of trw data.
recon_start_ppe <- min(df_list_ppe$start_years)
yearly_fig_ppe_df <- get_yearly_pred(fit_ppe_linear, clim_df_ppe, recon_start_ppe)

# get predictions and instrumental climate in weighted anomaly space
yearly_fig_ppe_df <- yearly_fig_ppe_df %>%
  mutate(
    recon_start = recon_start_ppe,
    calib_start = calib_start_ppe
  )

# plot reconstruction
fig_recon(yearly_fig_ppe_df,
          start_year_obs = 1600, calib_start_ppe,
          pred_name = "clim_mean_pred",
          truth_name = "mean_w_trueclim"
)[[2]] +
  geom_ribbon(aes(x = year, ymin = q5_w_clim, ymax = q95_w_clim), alpha = 0.5) +
  labs(title = str_glue('linear ppe using log(precip + 0.001, snr = {snr_ppe}'),
       subtitle = 'anomaly space',
       y = 'w_precip log anomaly (sds)'
       )


## repeat, but get the reconstruction in terms of DJF precipitation (original units mm/month)

djf_df_ppe <- process_climate(sim_precip_df,
                              clim_name = 'precip',
                              w_agg = w_djf, anom = F)

djf_fig_ppe_df <- get_djf_pred(fit_ppe_linear, djf_df_ppe, recon_start_ppe)
fig_recon(djf_fig_ppe_df,
          1600, calib_start_ppe,
          display_years = c(1750, 1950)
)[[2]] +
  labs(y = 'DJF mean precip (mm/month)')


